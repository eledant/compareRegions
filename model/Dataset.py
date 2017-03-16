from Data import Data
import re, random, pprint, copy

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
				self.order.append( region_name )

# -----------------------------------------------------------------------------------------

	####################################################
	### Calculate the regions position of the genome ###
	####################################################
	def calcGenomeCoords(self):
		pos = 0
		for region in self.order:
			data = self[region][0]
			data['startCoord'] = pos
			pos += data['chromEnd'] - data['chromStart'] + 1
			data['endCoord'] = pos - 1
	

	##################################################################
	### Calculate the regions coordinates of the 'A' and 'B' files ###
	##################################################################
	def calcDataCoords(self, genome):
		for region in self:		
			for data in self[region]:
				genomeData = genome[region][0]
				# Do something? Not in the first row?
				if 'strand' in genomeData :
					if genomeData['strand'] == '-':
						data['startCoord'] = genomeData['endCoord'] - (data['chromEnd'] - genomeData['chromStart'])
						if 'strand' in data:
							if data['strand'] == '+':
								data['coordStrand'] = '-'
							else:
								data['coordStrand'] = '+'
				else:
					data['startCoord'] = genomeData['startCoord'] + (data['chromStart'] - genomeData['chromStart'])
					if 'strand' in data:
						data['coordStrand'] = data['strand']

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
				randFile = copy.deepcopy(self)
				randFile.remapData(refGenome, randGenome)
				randFile.printRemapData(randGenome)
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

		# Restructurate the genome
		for chrom in self:
			if len(self[chrom]) == 2:
				if self[chrom][1]['chromEnd'] > self[chrom][0]['chromEnd']:
					self[chrom][0]['chromEnd'] = self[chrom][1]['chromEnd']
				else:
					self[chrom][0]['chromStart'] = self[chrom][1]['chromStart']
				del self[chrom][1]
		# Initialize variables
		lastCoord = -1
		regions = self.keys()
		random.shuffle(regions)
		firstRegion = regions.pop(0)
		data = self[firstRegion][0]
		data['strand'] = random.choice(['-', '+'])
		randomStart = random.randint( data['chromStart'], data['chromEnd']-1 )
		self[firstRegion].append( data.copy() )
		# Create random starting point
		if data['strand'] == '+':
			data['chromStart'] = randomStart
			self[firstRegion][1]['chromEnd'] = randomStart - 1
		else:
			data['chromEnd'] = randomStart 
			self[firstRegion][1]['chromStart'] = randomStart + 1
		lastCoord = remapGenome( data, lastCoord )		
		# Calculate new coordinates
		for region in regions:
			self[region][0]['strand'] = random.choice(['-', '+'])
			lastCoord = remapGenome( self[region][0], lastCoord )
		remapGenome( self[firstRegion][1], lastCoord )
		

	#####################################################
	### Randomize a file dataset with new coordinates ###
	#####################################################
	def remapData(self, refGenome, randGenome):
		for chrom in self:
			for data in self[chrom]:
				# Initialize variables
				refStartCoord = refGenome[chrom][0]['startCoord']
				randDataGenome = randGenome[chrom][0]
				diffStart = 0
				data['oldStartCoord'] = data['startCoord']
				# If data overlap the split chromosome
				if (data['chromStart'] < randDataGenome['chromEnd'] and data['chromEnd'] > randDataGenome['chromEnd']) or (data['chromStart'] < randDataGenome['chromStart'] and data['chromEnd'] > randDataGenome['chromStart']):
					data['tag'] = 'split'
					if randDataGenome['strand'] == '+':
						data['startCoord'] = randGenome[chrom][1]['startCoord'] + (data['startCoord'] - refStartCoord)
						data['endCoord'] = randGenome[chrom][0]['endCoord'] - (refGenome[chrom][0]['chromEnd'] - data['chromEnd'])
					else:
						data['startCoord'] = randGenome[chrom][0]['startCoord'] + (data['startCoord'] - refStartCoord)
						data['endCoord'] = randGenome[chrom][1]['endCoord'] - (refGenome[chrom][0]['chromEnd'] - data['chromEnd'])
				else:
					# If the chromosome is split in two
					if len(randGenome[chrom]) == 2:
						# If the data are on the other part of the split chromosome	
						if data['chromStart'] < randDataGenome['chromStart'] or data['chromEnd'] > randDataGenome['chromEnd']:
							randDataGenome = randGenome[chrom][1]
						diffStart = randDataGenome['chromStart'] - refGenome[chrom][0]['chromStart']
					# Calculate new coordinates
					data['startCoord'] = randDataGenome['startCoord'] + (data['startCoord'] - refStartCoord) - diffStart
					data['endCoord'] = data['startCoord'] + (data['chromEnd'] - data['chromStart'] )
				# Swap strand if data and the genome don't have the same 	### CORRECT????
				if randDataGenome['strand'] != data['strand']:
					data['coordStrand'] = randDataGenome['strand']

# -----------------------------------------------------------------------------------------

	##################################
	### Print the randomize genome ###
	##################################
	def printRandomizeGenome(self):
		print "Name\tchromStart\tchromEnd\tstrand\tstartCoord\tendCoord"
		print "########################################################################"
		for chrom in self:
			for genome in self[chrom]:
				print chrom, "\t", genome['chromStart'], "\t", genome['chromEnd'], "\t", genome['strand'], "\t", genome['startCoord'], "\t", genome['endCoord']
		print "########################################################################"


	##################################
	### Print the remapData output ###
	##################################
	def printRemapData(self, randGenome):
		print "Name\tchromStart\tchromEnd\tstrand\tstartCoord\tendCoord\tcoordStrand\toldStartCoord"
		print "##########################################################################"
		for chrom in self:
			for randData in randGenome[chrom]:
				print chrom, "\t", randData['chromStart'], "\t", randData['chromEnd'], "\t", randData['strand'], "\t", randData['startCoord'], "\t", randData['endCoord']
			for data in self[chrom]:
				print data['name'], "\t", data['chromStart'], "\t", data['chromEnd'], "\t", data['strand'], "\t", data['startCoord'], "\t", data['endCoord'], "\t", data['coordStrand'], "\t", data['oldStartCoord']
			print "##########################################################################"

