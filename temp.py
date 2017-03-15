				coordStrand = '+'
				if 'coordStrand' in data:
					coordStrand = data['coordStrand']

done = False
				for chromGenome in randGenome:		
					for dataGenome in randGenome[chromGenome]:
						if dataGenome['endCoord'] > startCoord and not done:
							newChrom = chromGenome
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
							newData[newChrom][-1]['oldChrom'] = chrom
							done = True
							
							#if endCoord > dataGenome['endCoord'] and dataGenome['strand']:
							#	if dataGenome['strand'] == '-':
							#	 	newData[newChrom][-1]['chromStart'] = dataGenome['chromStart']
							#	else:
							#		newData[newChrom][-1]['chromEnd'] = dataGenome['chromEnd']
							#	startCoord = (dataGenome['endCoord'] + 1) % self.size()
							#	overhang = endCoord - dataGenome['endCoord']
							#	endCoord = startCoord + overhang - 1
