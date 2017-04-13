from Data import Data
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
				self[region_name].append( Data(values, attributes) )
				if region_name not in self.order:
					self.order.append( region_name )
		for chrom in self.order:
			self[chrom] = sorted(self[chrom], key=lambda k: k['chromStart']) 

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
		# By default : Shuffle dataset regions et calculate new coordinates according to a new starting point
		if args['-l'] == 'def':
			randGenome.randomizeGenome()				
		# Calculate new coordinates of the data according to the randomize genome
		randFile = copy.deepcopy(self)
		randFile.remapData(randGenome)
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
			# Replace the lastChrom in a randomly order
			randPos = random.randint( 1, len(chroms)-1 )
			self.order = self.order[:randPos] + [lastChrom] + self.order[randPos:]
			chroms = chroms[:randPos] + [lastChrom] + chroms[randPos:]

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
		if 'split' in randGenome.order[-1]:
			self.order += [ randGenome.order[-1] ]
		# Initialize variables
		firstChrom = self.order[0]
		lastChrom = self.order[-1]
		toDelete = []
		# Add the split chromosome to the dataset
		self[lastChrom]= []
		for chrom in self.order[:-1]:
			for data in self[chrom]:
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
		overA, overB, over_bp, i, startIndex, scoreRegionB, nameRegionA, nameRegionB  = 0, 0, 0, 0, 0, [], [], []
		regionsA = getRegions(self)
		regionsB = getRegions(datasetB)
		# Loop on all the data in datasetA
		for dataA in regionsA:
			startA = dataA['startCoord']				
			endA = dataA['endCoord']
			scoreA = 1
			if args['-i'] not in ['A', 'AB'] and 'score' in dataA:
				scoreA = dataA['score']
			i = startIndex
			indexPossible = True
			# Loop on data in datasetB according to the index
			while i < len(regionsB):
				dataB = regionsB[i]
				startB = dataB['startCoord']
				endB = dataB['endCoord']
				# If there are an overlap bewteen A and B
				if endA >= startB and endB >= startA:
					scoreB = 1
					if args['-i'] not in ['B', 'AB'] and 'score' in dataB:
						scoreB = dataB['score']
					# Calculate the number of overlap
					if scoreA * scoreB >= 1:
						if indexPossible:
							if 'split' in dataA:
								nameRegionA.append(dataA['name'])
							else:
								overA += 1
						if 'split' in dataB:
							nameRegionB.append(dataB['name'])
						else:
							scoreRegionB.append(i)
					# Calculate the size of the overlap
					over_bp += getOverlap(startA, endA, startB, endB) * scoreA * scoreB
					indexPossible = False
				# Set index if possible (condition + never get an overlap for this dataA)
				elif startA > endB and indexPossible:
					startIndex = i + 1
				# Break if it's impossible to find overlap
				elif endA < startB:
					break
				i += 1
		
		overA += len(list(set(nameRegionA)))
		overB = len(list(set(scoreRegionB))) + len(list(set(nameRegionB)))
		return [overA, overB, over_bp]

# -----------------------------------------------------------------------------------------

	#####################################
	### Print a entire Dataset object ###
	#####################################
	def printDataset(self, filename):
		path = 'test_Files/%s' % filename
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
					score = "\t%.2f" % data['score']
				FH.write("%s\t%d\t%d%s%s%s\n" % (chrom, data['startCoord'], data['endCoord'], name, score, strand) )
		FH.close()



