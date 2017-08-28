#!/usr/bin/python
# -*- coding: utf-8 -*-
import random

# Parameters
folderPath = './'
name = 'datasetB_'+str(regionMax)+'.bed'
lines = []
nbRegionList = ['1', '10', '25', '50', '75', '100', '125', '150', '175', '200', '225', '250', '500', '750', '1000', '1250', '1500', '1750', '2000']

# Open reference dataset to sample
F = open(folderPath+name, 'r')
for line in F:
	lines.append(line)
# Delete header
lines.pop(0)
F.close()

regionMax= len(lines)

# For each dataset created
for nbRegion in nbRegionList:
	# Generate sampled dataset
	alreadyDone = []
	name = 'datasetB_'+str(nbRegion)+'.bed'
	F = open(folderPath+name, 'w')
	F.write('track\ttype=Bed\tname="%s"\n' % name)
	for indexLine in range(int(nbRegion)):
		while True:
			randRegion = random.randint(0, regionMax-1)
			if randRegion not in alreadyDone:
				F.write(lines[randRegion])
				alreadyDone.append(randRegion)
				break
	F.close()


