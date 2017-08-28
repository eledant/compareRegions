#!/usr/bin/python
# -*- coding: utf-8 -*-
import random, re

importPath = './'
importFile = 'datasetB'
genomePath = '../'
genomeName = 'genome_test.bed'
dictRegion = {}
xTime = 10


# Save regions of an import fileA
F = open(importPath+importFile+'.bed', 'r')
for line in F:
	if re.search('^track', line) is None:
		line = line.rstrip('\n')
		chrmName, start, end, idRegion = re.split('\t', line)
		start, end, idRegion = int(start), int(end), idRegion
		if chrmName not in dictRegion:
			dictRegion[chrmName] = [ [start, end, idRegion] ]
		else:
			dictRegion[chrmName].append( [start, end, idRegion] )
	else:
		firstLine = line
F.close()

# Create new fileA
for x in range(1,xTime+1):
	F.write(firstLine)
	for chrmName in dictRegion:
		for i in range(len(dictRegion[chrmName])):
			start, end, idRegion = dictRegion[chrmName][i]
			regionStarts = []
			for nbRegion in range(x):
				regionStart = random.randint(start, end)
				regionStarts.append( regionStart )
			regionStarts = sorted(regionStarts)
			F.write('%s\t%d\t%d\t%s\n' %(chrmName, start, regionStarts[0]-1, idRegion+'_0'))
			for j in range(0,len(regionStarts)-1):
				F.write('%s\t%d\t%d\t%s\n' %(chrmName, regionStarts[j], regionStarts[j+1]-1, idRegion+'_'+str(j)))
			F.write('%s\t%d\t%d\t%s\n' %(chrmName, regionStarts[-1], end, idRegion+'_last'))
	F.close()
