# A dictionary with the attributes as key and the corresponding column as value 
class Data(dict):

	# Set the dictionary
	def __init__(self, values, attributes):
		length = min( len(values), len(attributes) )
		for i in range(length):
			if values[i].isdigit():
				self[ attributes[i] ] = int(values[i])
			else:
				self[ attributes[i] ] = values[i]
	
	# Return the size of the region
	def size(self):
		return self['chromEnd'] - self['chromStart']
