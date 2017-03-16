# A dictionary with the attributes as key and the corresponding column as value 
class Data(dict):

	# Set the dictionary
	def __init__(self, values, attributes):
		length = min( len(values), len(attributes) )
		for i in range(length):
			if str(values[i]).isdigit():
				self[ attributes[i] ] = int(values[i])
			else:
				self[ attributes[i] ] = values[i]

	# Change the coordinates
	def remapGenome(self, lastCoord):
		self['startCoord'] = lastCoord + 1
		self['endCoord'] = self['startCoord'] + (self['chromEnd'] - self['chromStart'])			
		return self['endCoord']
