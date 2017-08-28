#!/usr/bin/python
# -*- coding: utf-8 -*-
import random, re

importPath = './'
importFile = 'datasetB'
genomePath = '../../'
genomeName = 'genome_test.bed'
dictRegion = {}
nbRegions = 500
averageSize = 0
nbElement = 0
xTime = 10

# Save chroms coordinates and name of the genome file in a list
chrms = {}
coordRemap = 0
try:
	FG = open(genomePath+genomeName, 'r')
	for line in FG:
		if re.search('^track', line) is None:
			line = line.rstrip('\n')
			name, start, end = re.split('\t', line)
			chrms[name] = [int(start), int(end), coordRemap]
			coordRemap += int(end) - int(start)
except IOError:
    	print "Could not read file:", genomeName
    	exit()
finally:
   	FG.close()

# Save regions of an import file
F = open(importPath+importFile+'.bed', 'r')
for line in F:
	if re.search('^track', line) is None:
		line = line.rstrip('\n')
		chrmName, start, end, idRegion = re.split('\t', line)
		start, end, idRegion = int(start), int(end), idRegion
		averageSize += end - start + 1
		nbElement += 1 
		if chrmName not in dictRegion:
			dictRegion[chrmName] = [ [start, end, idRegion] ]
		else:
			dictRegion[chrmName].append( [start, end, idRegion] )
	else:
		firstLine = line
F.close()

averageSize /= nbElement

# Create new fileB
for x in range(1,xTime+1):
	F = open(importPath+importFile+'_'+str(x)+'.bed', 'w')
	F.write(firstLine)
	for chrmName in dictRegion:
		for i in range(len(dictRegion[chrmName])):
			start, end, idRegion = dictRegion[chrmName][i]
			F.write('%s\t%d\t%d\t%s\n' %(chrmName, start, end, idRegion))
			dictRegion[chrmName][i] = [start, end, idRegion]
	for idR in range(nbRegions):
		randomStart = random.randint(0, coordRemap)
		for chrmName in chrms:
			chrmStart, chrmEnd, startCoord = chrms[chrmName]
			if randomStart > startCoord and randomStart < (startCoord+chrmEnd-chrmStart):
				break
		randomSize = averageSize
		start = random.randint(chrms[chrmName][0], chrms[chrmName][1]-randomSize)
		end = start + randomSize - 1
		idRegion = 'noise_%d_%d' %(x, idR)
		F.write('%s\t%d\t%d\t%s\n' %(chrmName, start, end, idRegion))
		if chrmName not in dictRegion:
			dictRegion[chrmName] = [ [start, end, idRegion] ]
		else:
			dictRegion[chrmName].append([start, end, idRegion])
	F.close()
