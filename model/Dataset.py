from OverlapMatrix import OverlapMatrix
import re, random, pprint, copy, datetime, os.path

# A dictonary with the chromosom number as key and a list of Data as value
class Dataset(dict):

	#####################################################
	### Open the file and call the 'getData' function ###
	#####################################################
	def __init__(self, filename=None):
		if filename:
			self.name = filename
			self.order = []
			try:
				FILE = open(filename, 'r') 
				self.getData(FILE)
			except IOError:
    				print "Could not read file:", filename
    				exit()
			finally:
   				FILE.close()

	
	################################################################
	### Extract Data from 'FILE', sort them by the 'region_name' ###
	################################################################
	def getData(self, FILE):
		# Set the dictionary
		def setDict(values, attributes):
			dic = {}
			length = min( len(values), len(attributes) )
			for i in range(length):
				# Special case when option -l rdg is set (TO CHANGE for better maintenance)
				if attributes[i] == 'score':
					dic[ attributes[i] ] = float(values[i])
				elif str(values[i]).isdigit():
					dic[ attributes[i] ] = int(values[i])
				else:
					dic[ attributes[i] ] = values[i]
			return dic

		attributes = ['chromStart', 'chromEnd']
		# Check if it's a BED or BEDGRAPH file
		for line in FILE:
			if re.search('^track', line) is not None:
				if re.search('bedgraph', line, re.I) is not None:
					attributes.append('score')
				else:
					attributes.extend( ['name', 'score', 'strand'] )
				break
		# Get the values in each column
		for line in FILE:
			if re.search('^#', line) is None:
				line = line.rstrip('\n')
				values = re.split('\t', line)
				region_name = values.pop(0)
				if region_name not in self:
					self[region_name] = []
				self[region_name].append( setDict(values, attributes) )
				if region_name not in self.order:
					self.order.append( region_name )
		# Sort regions by order of height
		for chrom in self.order:
			self[chrom] = sorted(self[chrom], key=lambda k: k['chromStart'])


	#################################################
	### Return the size and the number of regions ###
	#################################################
	def getSize(self):
		nbRegions, size = 0, 0
		for chrom in self:
			for data in self[chrom]:
				nbRegions += 1
				size += data['chromEnd'] - data['chromStart'] + 1
		return [nbRegions, size]
		
# -----------------------------------------------------------------------------------------

	####################################################
	### Calculate the regions position of the genome ###
	####################################################
	def calcGenomeCoords(self):
		pos = 0
		for chrom in self.order:
			data = self[chrom][0]
			data['startCoord'] = pos
			pos += data['chromEnd'] - data['chromStart'] + 1
			data['endCoord'] = pos - 1
	

	##################################################################
	### Calculate the regions coordinates of the 'A' and 'B' files ###
	##################################################################
	def calcDataCoords(self, genome):
		for chrom in self:		
			for data in self[chrom]:
				genomeData = genome[chrom][0]
				data['startCoord'] = genomeData['startCoord'] + (data['chromStart'] - genomeData['chromStart'])
				data['endCoord'] = data['startCoord'] + (data['chromEnd'] - data['chromStart'])

# -----------------------------------------------------------------------------------------

	##############################################
	### Create n randomizations of the dataset ###
	##############################################
	def randomize(self, refGenome, args, filename):
		randGenome = copy.deepcopy(refGenome)
		if args['-r']:
			random.seed()
		randFile = copy.deepcopy(self)
		# 1=Standard ; 2=Within chroms ; 3=Within rand
		randType = 1
		# By default : Shuffle dataset regions et calculate new coordinates according to a new starting point
		if randType == 1:
			randGenome.randomizeGenome()
			# Calculate new coordinates of the data according to the randomize genome
			randFile.remapData(randGenome)
		elif randType == 2:
			randGenome.randomizeGenomeWithinChrom()
			randFile.remapDataWithinChrom(randGenome)
		else:
			randGenome.randomizeGenomeWithinRandomly()
			randFile.remapDataRandomly(randGenome)		
		
		if args['-v'] in ['all', 'randG']: randGenome.printDataset(filename)
		return randFile


	##################################
	### Randomize a genome dataset ###
	##################################
	def randomizeGenome(self):
		# Change the coordinates
		def remapGenome(data, lastCoord):
			data['startCoord'] = lastCoord + 1
			data['endCoord'] = data['startCoord'] + (data['chromEnd'] - data['chromStart'])
			return data['endCoord']

		# Initialize variables
		lastCoord = -1
		chroms = self.order
		# Create random starting point
		randomStart = random.randint( 0, self[chroms[-1]][0]['endCoord'] )
		random.shuffle(chroms)
		for i in range(len(chroms)):
			data = self[chroms[i]][0]
			if data['startCoord'] <= randomStart <= data['endCoord']:
				chroms = chroms[i:] + chroms[:i]
				break

		firstChrom = chroms.pop(0)
		data = self[firstChrom][0]
		data['strand'] = random.choice(['-', '+'])
		randomStart = data['chromStart'] + ( randomStart - data['startCoord'] )
		# Save the order as index for later comparaison
		self.order = [firstChrom] + chroms[:]
		# Add the second part of the split chrom
		if randomStart != data['chromStart']:
			lastChrom = firstChrom + '_split'
			self[lastChrom] = [data.copy()]
			self[lastChrom][0]['chromEnd'] = randomStart - 1
			self.order.append(lastChrom)
			chroms.append(lastChrom)
		data['chromStart'] = randomStart
		lastCoord = remapGenome( data, lastCoord )
		# Calculate new coordinates
		for chrom in chroms:
			self[chrom][0]['strand'] = random.choice(['-', '+'])
			lastCoord = remapGenome( self[chrom][0], lastCoord )
		

	#####################################################
	### Randomize a file dataset with new coordinates ###
	#####################################################
	def remapData(self, randGenome):
		# Set the new chromosome order in the randomize dataset
		temp_chrom = []
		for chrom in randGenome.order:
			if chrom in self.order:
				temp_chrom.append(chrom)
		self.order = temp_chrom[:]
		# Add the split chromosome to the dataset
		if 'split' in randGenome.order[-1]:
			lastChrom = randGenome.order[-1]
			self[lastChrom]= []
		# Initialize variables
		firstChrom = self.order[0]
		toDelete = []
		for chrom in self.order:
			for i in range(len(self[chrom])):
				data = self[chrom][i]
				randDataGenome = randGenome[chrom][0]
				# If the chromosome is split in two
				if chrom == firstChrom and data['chromStart'] < randDataGenome['chromStart']:
					data_split = data.copy()
					# If data overlap the split chromosome
					if data['chromEnd'] >= randDataGenome['chromStart']:
						# Modify the data part on the firstChrom
						data['chromStart'] = randDataGenome['chromStart']
						data['startCoord'] = 0
						data['endCoord'] = data['startCoord'] + (data['chromEnd'] - data['chromStart'] )
						data['split'] = True
						# Modify the data part on the lastChrom
						data_split['chromEnd'] = randGenome[lastChrom][0]['chromEnd']
						data_split['startCoord'] = randGenome[lastChrom][0]['endCoord'] - (data_split['chromEnd'] - data_split['chromStart'])
						data_split['endCoord'] = randGenome[lastChrom][0]['endCoord']
						data_split['split'] = True
					else:
						# Modify the data which is only on the lastChrom
						data_split['startCoord'] = randGenome[lastChrom][0]['startCoord'] + (data_split['chromStart'] - randGenome[lastChrom][0]['chromStart'])
						data_split['endCoord'] = data_split['startCoord'] + (data_split['chromEnd'] - data_split['chromStart'] )
						# Save the firstChrom data to delete it later (don't exist now)
						toDelete.append( data )
					self[lastChrom].append(data_split)	
				else:
					# Calculate new coordinates
					data['startCoord'] = randDataGenome['startCoord'] + (data['chromStart'] - randDataGenome['chromStart'])
					data['endCoord'] = data['startCoord'] + (data['chromEnd'] - data['chromStart'] )
		# Add the split chromosome to the order list
		if 'split' in randGenome.order[-1]:
			self.order.append( lastChrom )
		# Delete the useless regions in the firstChrom
		self[firstChrom] = [data for data in self[firstChrom] if data not in toDelete]


		# Inverse coordinates if '-' strand
		#for chrom in self.order:
		#	randDataGenome = randGenome[chrom][0]
		#	if randDataGenome['strand'] == '-':
		#		for data in self[chrom]:
		#			data['endCoord'] = randDataGenome['endCoord'] - (data['startCoord'] - randDataGenome['startCoord'])
		#			data['startCoord'] = data['endCoord'] - (data['chromEnd'] - data['chromStart'])
	

# -----------------------------------------------------------------------------------------
	
	########################################################################
	### Count overlaps between regions of this A dataset and a B dataset ###
	########################################################################
	def compareData(self, datasetB, args):
		# Return the overlap size		
		def getOverlap(start1, end1, start2, end2):
			return min(end1, end2) - max(start1, start2) + 1
		# Return a list of all the data
		def getRegions(dataset):
			regions = []
			for chrom in dataset:
				regions += dataset[chrom]
			regions = sorted(regions, key=lambda k: k['startCoord']) 
			return regions
			
		# Initialize variables
		pointerB, splitRegionA, splitRegionB  = 0, [], []
		regionsA, regionsB = getRegions(self), getRegions(datasetB)
		nA, sizeA = self.getSize()
		nB, sizeB = datasetB.getSize()
		matrix = OverlapMatrix(nA, nB, sizeA, sizeB)
		matrix.getPearsonVariables(regionsA, regionsB)
		# Loop on all the data in datasetA
		for indexA in range(len(regionsA)):
			dataA = regionsA[indexA]
			startA, endA, scoreA = dataA['startCoord'], dataA['endCoord'], 1
			if args['-i'] not in ['A', 'AB'] and 'score' in dataA:
				scoreA = dataA['score']
			indexB = pointerB
			indexPossible = True
			# Loop on data in datasetB according to the index
			while indexB < len(regionsB):
				dataB = regionsB[indexB]
				startB, endB = dataB['startCoord'], dataB['endCoord']
				# If there are an overlap bewteen A and B
				if endA >= startB and endB >= startA:
					scoreB = 1
					if args['-i'] not in ['B', 'AB'] and 'score' in dataB:
						scoreB = dataB['score']
					# Save split regions names
					if 'split' in dataA:
						splitRegionA.append(dataA['name'])
					if 'split' in dataB:
						splitRegionB.append(dataB['name'])
					# Calculate the size of the overlap
					overlap = getOverlap(startA, endA, startB, endB) * scoreA * scoreB
					# Save the results in the matrix
					matrix.addOverlap( overlap, indexA, indexB )
					indexPossible = False
				# Set index if possible (condition + never get an overlap for this dataA)
				if startA > endB and indexPossible:
					pointerB = indexB + 1
				# Break if it's impossible to find overlap
				elif endA < startB:
					break
				indexB += 1

		# Add split regions to the matrix object
		matrix.addSplitRegions( splitRegionA, splitRegionB )
		matrix.setOverlapMeasures()
		return matrix

# -----------------------------------------------------------------------------------------

	#####################################
	### Print a entire Dataset object ###
	#####################################
	def printDataset(self, filename):
		path = 'testFiles/%s.bed' % filename
		FH = open(path, 'w')
		BEDname = re.split('/', self.name)
		FH.write('track\ttype=Bed\tname="%s"\n' % BEDname[-1])
		for chrom in self.order:
			for data in self[chrom]:
				name, strand, score = '', '', ''
				if 'name' in data:
					name = "\t%s" % data['name']
				if 'strand' in data:
					strand = "\t%s" % data['strand']
				if 'score' in data:
					score = "\t%g" % data['score']
				FH.write("%s\t%d\t%d%s%s%s\n" % (chrom, data['startCoord'], data['endCoord'], name, score, strand) )
		FH.close()


# -----------------------------------------------------------------------------------------

	##################################
	### Randomize a genome dataset ###
	##################################
	def randomizeGenomeWithinChrom(self):
		# Change the coordinates
		def remapGenome(data, lastCoord):
			data['startCoord'] = lastCoord + 1
			data['endCoord'] = data['startCoord'] + (data['chromEnd'] - data['chromStart'])
			return data['endCoord']

		newOrder = []
		lastCoord = -1
		# Create random starting point for each chrom
		for i in range(len(self.order)):
			# Initialize variables
			chrom = self.order[i]
			data = self[chrom][0]
			# Create random starting point
			randomStart = random.randint( data['chromStart'], data['chromEnd'] )
			data['strand'] = random.choice(['-', '+'])
			# Create new chrom with new coordinates
			newChrom = chrom + '_split'
			self[newChrom] = [data.copy()]
			self[newChrom][0]['chromEnd'] = randomStart - 1
			# Save the order as index for later comparaison
			newOrder.append(chrom)
			newOrder.append(newChrom)
			# Change old chrom coordinates
			data['chromStart'] = randomStart
			# Remap chrom coordinates
			lastCoord = remapGenome( self[chrom][0], lastCoord )
			lastCoord = remapGenome( self[newChrom][0], lastCoord )
		# Add new chroms
		self.order = newOrder[:]


	#####################################################
	### Randomize a file dataset with new coordinates ###
	#####################################################
	def remapDataWithinChrom(self, randGenome):
		# Set the new chromosome order in the randomize dataset
		temp_chrom, i = [], 0
		while i < len(randGenome.order):
			if randGenome.order[i] in self.order:
				temp_chrom.append( randGenome.order[i] )
				temp_chrom.append( randGenome.order[i+1] )
			i += 2
		self.order = temp_chrom[:]

		# Initialize variables
		i = 0
		while i < len(self.order):
			chrom = self.order[i]
			chromSplit = self.order[i+1]
			self[chromSplit] = []
			toDelete = []
			for j in range(len(self[chrom])):
				data = self[chrom][j]
				randDataGenome = randGenome[chrom][0]
				# If data overlap the split chromosome
				if data['chromEnd'] >= randDataGenome['chromStart'] and data['chromStart'] < randDataGenome['chromStart']:
					data_split = data.copy()
					# Modify the data part on the firstChrom
					data['chromStart'] = randDataGenome['chromStart']
					data['startCoord'] = randDataGenome['startCoord']
					data['endCoord'] = data['startCoord'] + (data['chromEnd'] - data['chromStart'] )
					data['split'] = True
					# Modify the data part on the splitChrom
					data_split['chromEnd'] = randGenome[chromSplit][0]['chromEnd']
					data_split['endCoord'] = randGenome[chromSplit][0]['endCoord']
					data_split['startCoord'] = data_split['endCoord'] - (data_split['chromEnd'] - data_split['chromStart'])
					
					data_split['split'] = True
					self[chromSplit].append(data_split)
				else:
					# If data on the split chrom
					if data['chromEnd'] < randDataGenome['chromStart']:
						# Modify the data which is only on the splitChrom
						data['startCoord'] = randGenome[chromSplit][0]['startCoord'] + (data['chromStart'] - randGenome[chromSplit][0]['chromStart'])
						data['endCoord'] = data['startCoord'] + (data['chromEnd'] - data['chromStart'] )
						self[chromSplit].append(data)
						toDelete.append(j)
					else:
						# Calculate new coordinates
						data['startCoord'] = randDataGenome['startCoord'] + (data['chromStart'] - randDataGenome['chromStart'])
						data['endCoord'] = data['startCoord'] + (data['chromEnd'] - data['chromStart'] )		
			i += 2
			# Delete the useless regions in the chrom
			x = 0
			for e in toDelete:
				del self[chrom][e-x]	
				x += 1

# -----------------------------------------------------------------------------------------

	##################################
	### Randomize a genome dataset ###
	##################################
	def randomizeGenomeWithinRandomly(self):
		# Change the coordinates
		def remapGenome(data, lastCoord):
			data['startCoord'] = lastCoord + 1
			data['endCoord'] = data['startCoord'] + (data['chromEnd'] - data['chromStart'])
			return data['endCoord']

		lastCoord = -1
		# Create random starting point for each chrom
		for i in range(len(self.order)):
			# Initialize variables
			chrom = self.order[i]
			# Remap chrom coordinates
			lastCoord = remapGenome( self[chrom][0], lastCoord )


	#####################################################
	### Randomize a file dataset with new coordinates ###
	#####################################################
	def remapDataRandomly(self, randGenome):
		# Set the new chromosome order in the randomize dataset
		temp_chrom, i = [], 0
		while i < len(randGenome.order):
			if randGenome.order[i] in self.order:
				temp_chrom.append( randGenome.order[i] )
			i += 1
		self.order = temp_chrom[:]

		# Initialize variables
		i = 0
		while i < len(self.order):
			chrom = self.order[i]
			for j in range(len(self[chrom])):
				data = self[chrom][j]
				randDataGenome = randGenome[chrom][0]
				# Calculate new coordinates
				data['startCoord'] = random.randint( randDataGenome['startCoord'], randDataGenome['endCoord'] - (data['chromEnd'] - data['chromStart']) )
				data['endCoord'] = data['startCoord'] + (data['chromEnd'] - data['chromStart'])		
			i += 1


