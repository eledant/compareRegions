#!/usr/bin/python
# -*- coding: utf-8 -*-
import random, re

option = '_0'
importPath = './'
importFile = 'datasetA'
genomePath = '../../'
genomeName = 'genome_test.bed'
chrms, dictRegion = {}, {}

# Save chroms coordinates and name of the genome file in a list
try:
	FG = open(genomePath+genomeName, 'r')
	for line in FG:
		if re.search('^track', line) is None:
			line = line.rstrip('\n')
			name, start, end = re.split('\t', line)
			chrms[name] = [int(start), int(end)]
except IOError:
    	print "Could not read file:", genomeName
    	exit()
finally:
   	FG.close()


# Save regions of an import fileA
F = open(importPath+importFile+option+'.bed', 'r')
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


# Sort regions
for chrmName in dictRegion:
	dictRegion[chrmName].sort(key=lambda x: x[0])


# Create new fileB
F = open('datasetB'+option+'.bed', 'w')
F.write(firstLine)
for chrmName in chrms:
	startChrm, endChrm = chrms[chrmName]
	if chrmName in dictRegion:
		for i in range(len(dictRegion[chrmName])):
			startRegion, endRegion, idRegion = dictRegion[chrmName][i]
			idRegion = "DataB_%s_%d" %(chrmName, i)
			if i == 0:
				F.write('%s\t%d\t%d\t%s\n' %(chrmName, startChrm, startRegion-1, idRegion))
			else:
				F.write('%s\t%d\t%d\t%s\n' %(chrmName, previousEnd+1, startRegion-1, idRegion))
			previousEnd = endRegion
		F.write('%s\t%d\t%d\t%s\n' %(chrmName, previousEnd+1, endChrm, idRegion))
	else:
		idRegion = "DataB_%s" %chrmName
		F.write('%s\t%d\t%d\t%s\n' %(chrmName, startChrm, endChrm, idRegion))
F.close()



