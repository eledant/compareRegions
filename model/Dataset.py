from Data import Data
import re, random, pprint, copy

# A dictonary with the chromosom number as key and a list of Data as value
class Dataset(dict):

	# Open the file and call the 'getData' function
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

	
	
	# Extract Data from 'FILE', sort them by the 'region_name'
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
				

	# Calculate the regions position of the genome
	def calcGenomeCoords(self):
		pos = 0
		for region in self.order:
			data = self[region][0]
			data['startCoord'] = pos
			pos += data['chromEnd'] - data['chromStart'] + 1
			data['endCoord'] = pos - 1
	

	# Calculate the regions coordinates of the 'A' and 'B' files
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


	# Create n randomizations of the dataset
	def randomize(self, args, refGenome):
		randGenome = copy.deepcopy(refGenome)
		if args['-r']:
			random.seed()
		for nbRandom in range( int(args['-n']) ):
			if args['-m'] == 'rel':
				# Copy
				randGenome.randomizeGenome()
				self.remapData(randGenome)


	# Randomize a genome dataset
	def randomizeGenome(self):
		# Change the coordinates
		def remapGenome(data, lastCoord):
			data['startCoord'] = lastCoord + 1
			data['endCoord'] = data['startCoord'] + (data['chromEnd'] - data['chromStart'])			
			return data['endCoord']

		lastCoord = -1
		regions = self.keys()
		random.shuffle(regions)
		firstRegion = regions.pop(0)

		data = self[firstRegion][0]
		data['strand'] = random.choice(['-', '+'])
		randomStart = random.randint( data['chromStart'], data['chromEnd'] )
		self[firstRegion].append( data.copy() )
		if data['strand'] == '+':
			data['chromStart'] = randomStart
			self[firstRegion][1]['chromEnd'] = randomStart - 1
		else:
			data['chromEnd'] = randomStart
			self[firstRegion][1]['chromStart'] = randomStart + 1
		lastCoord = remapGenome( data, lastCoord )		

		for region in regions:
			self[region][0]['strand'] = random.choice(['-', '+'])
			lastCoord = remapGenome( self[region][0], lastCoord )
		remapGenome( self[firstRegion][1], lastCoord )
		

	# Randomize a file dataset with new coordinates
	def remapData(self, randGenome):
		newData = Dataset()
		for region in self:		
			for data in self[region]:
				startCoord = data['startCoord']
				endCoord = startCoord + ( data['chromEnd'] - data['chromStart'] )
				coordStrand = '+'
				if 'coordStrand' in data:
					coordStrand = data['coordStrand']
				
				done = False
				for region in randGenome:		
					for dataGenome in randGenome[region]:
						if dataGenome['endCoord'] > startCoord and not done:
							newChrom = region
							if 'strand' in dataGenome and dataGenome['strand'] == '-':
								chromEnd = dataGenome['chromEnd'] - (startCoord - dataGenome['startCoord'])
								chromStart = chromEnd - (endCoord - startCoord)
								if coordStrand == '+':
									strand = '-'
								else:
										strand = '+'
							else:
								chromStart = dataGenome['chromStart'] + (startCoord - dataGenome['startCoord'])
								chromEnd = chromStart + (endCoord - startCoord)
								strand = coordStrand
							if newChrom not in newData:
								newData[newChrom] = []	
							values = ['chromStart', 'chromEnd', 'score', 'strand', 'startCoord', 'coordStrand']
							attributes = [chromStart, chromEnd, data['score'], strand, startCoord, coordStrand ] 				
							newData[newChrom].append( Data(attributes, values) )
							newData[newChrom][-1]['oldChrom'] = region
							done = True
							
							#if endCoord > dataGenome['endCoord'] and dataGenome['strand']:
							#	if dataGenome['strand'] == '-':
							#	 	newData[newChrom][-1]['chromStart'] = dataGenome['chromStart']
							#	else:
							#		newData[newChrom][-1]['chromEnd'] = dataGenome['chromEnd']
							#	startCoord = (dataGenome['endCoord'] + 1) % self.size()
							#	overhang = endCoord - dataGenome['endCoord']
							#	endCoord = startCoord + overhang - 1
		pprint.pprint(newData)
			

	def size(self):
		size = 0
		for region in self:		
			for data in self[region]:
				size += data['chromEnd'] - data['chromStart'] +1



		



