from Data import Data
import re, random, pprint, copy, datetime

# A dictonary with the chromosom number as key and a list of Data as value
class Dataset(dict):

	#####################################################
	### Open the file and call the 'getData' function ###
	#####################################################
	def __init__(self, filename=None):
		if filename:
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
	def randomize(self, args, refGenome):
		randGenome = copy.deepcopy(refGenome)
		if args['-r']:
			random.seed()
		for nbRandom in range( int(args['-n']) ):
			# Preserve relational structure at all sub-region scales (default): permute dataset regions ### METHOD CHANGE
			if args['-m'] == 'rel':
				randGenome.randomizeGenome()
				if args['-v']: randGenome.printRandomizeGenome()
		randFile = copy.deepcopy(self)
		randFile.remapData(randGenome)
		if args['-v']: randFile.printRemapData(randGenome)
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
		# Restructure the genome if a chromosome is split in two
		def restructureGenome(self):
			firstChrom = self.order[0]
			lastChrom = self.order[-1]
			if '_split' in lastChrom:
				self[firstChrom][0]['chromStart'] = self[lastChrom][0]['chromStart']
				del self[lastChrom]
				del self.order[-1]
		
		restructureGenome(self)
		# Initialize variables
		lastCoord = -1
		chroms = self.keys()
		random.shuffle(chroms)
		firstChrom = chroms.pop(0)
		lastChrom = firstChrom + '_split'
		data = self[firstChrom][0]
		self[lastChrom] = [data.copy()]
		# Save the order as index for later comparaison
		self.order = [firstChrom] + chroms[:] + [lastChrom]
		# Create random starting point
		data['strand'] = random.choice(['-', '+'])
		randomStart = random.randint( data['chromStart'], data['chromEnd']-1 )
		data['chromStart'] = randomStart
		self[lastChrom][0]['chromEnd'] = randomStart - 1
		lastCoord = remapGenome( data, lastCoord )		
		# Calculate new coordinates
		chroms += [lastChrom]
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
			if chrom in self:
				for data in self[chrom]:
					randDataGenome = randGenome[chrom][0]
					# If the chromosome is split in two
					if chrom == firstChrom and data['chromStart'] < randDataGenome['chromStart']:
						data['split'] = True
						data_split = data.copy()
						# If data overlap the split chromosome
						if data['chromEnd'] >= randDataGenome['chromStart']:
							# Modify the data part on the firstChrom
							data['chromStart'] = randDataGenome['chromStart']
							data['startCoord'] = 0
							data['endCoord'] = data['startCoord'] + (data['chromEnd'] - data['chromStart'] )
							# Modify the data part on the lastChrom
							data_split['chromEnd'] = randGenome[lastChrom][0]['chromEnd']
							data_split['startCoord'] = randGenome[lastChrom][0]['endCoord'] - (data_split['chromEnd'] - data_split['chromStart'])
							data_split['endCoord'] = randGenome[lastChrom][0]['endCoord']
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
			for chrom in dataset.order:
				regions += dataset[chrom]
			regions = sorted(regions, key=lambda k: k['startCoord']) 
			return regions
		
		# Initialize variables	
		overA, overB, over_bp, i, startIndex, scoreRegionB  = 0, 0, 0, 0, 0, []
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
							overA += 1
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
		overB = len(list(set(scoreRegionB)))
		return [overA, overB, over_bp]


# -----------------------------------------------------------------------------------------

	##################################
	### Print the reference genome ###
	##################################
	def printReferenceGenome(self):
		print "###REFERENCE GENOME###"
		print "Name\tchromStart\tchromEnd\tstartCoord"
		print "##################################################"
		for chrom in self.order:
			for genome in self[chrom]:
				print chrom, "\t", genome['chromStart'], "\t", genome['chromEnd'], "\t",  genome['startCoord'], '\t', genome['name']
		print "##################################################"


	##################################
	### Print the randomize genome ###
	##################################
	def printRandomizeGenome(self):
		print "###RANDOMIZE GENOME###"
		print "Name\tchromStart\tchromEnd\tstrand\tstartCoord\tendCoord"
		print "########################################################################"
		for chrom in self.order:
			for genome in self[chrom]:
				print chrom, "\t", genome['chromStart'], "\t", genome['chromEnd'], "\t", genome['strand'], "\t", genome['startCoord'], "\t", genome['endCoord']
		print "########################################################################"


	##################################
	### Print the remapData output ###
	##################################
	def printRemapData(self, randGenome):
		print "###REMAP DATA###"
		print "Name\tchromStart\tchromEnd\tstrand\tstartCoord\tendCoord\tsplit"
		print "#####################################################################################"
		for chrom in self.order:
			for randData in randGenome[chrom]:
				print ">",chrom, "\t\t", randData['chromStart'], "\t", randData['chromEnd'], "\t", randData['strand'], "\t", randData['startCoord'], "\t", randData['endCoord']
			for data in self[chrom]:
				split = ''
				if 'split' in data :
					split = data['split']
				print data['name'], "\t", data['chromStart'], "\t", data['chromEnd'], "\t", data['strand'], "\t", data['startCoord'], "\t", data['endCoord'], "\t", split
			print "#####################################################################################"

