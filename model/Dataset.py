from Data import Data
import re, random

# A dictonary with the chromosom number as key and a list of Data as value
class Dataset(dict):

	# Open the file and call the 'getData' function
	def __init__(self, filename):
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
		if args['-r']:
			random.seed()
		for nbRandom in range( args['-n'] ):
			if args['-m'] == 'rel':
				refGenome.randomizeGenome()


	def randomizeGenome(self):
		strand = ['-', '+']
		randomStart = random.randint(0, self.calcSize())

		regions = self.keys()
		random.shuffle(regions)
		for region in regions:		
			data = self[region][0]
			data['strand'] = strand[ random.randint(0,1) ]
			
				
	def calcSize():
		size = 0
		for region in self:		
			for data in self[region]:
				size += data['chromEnd'] - data['chromStart'] + 1
		return size





		



