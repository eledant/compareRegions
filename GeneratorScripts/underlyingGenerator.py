#!/usr/bin/python
# -*- coding: utf-8 -*-
import random, re

importPath = './'
importFile = 'datasetB'
genomePath = '../'
genomeName = 'genome_test.bed'
dictRegion = {}
bpTotal = 100
xTime = 10

# Save chroms coordinates and name of the genome file in a list
chrms = {}
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

# Save regions of an import fileB
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

# Create new fileB
for x in range(1,xTime+1):
	F = open(importPath+importFile+'_'+str(x)+'.bed', 'w')
	F.write(firstLine)
	for chrmName in dictRegion:
		for i in range(len(dictRegion[chrmName])):
			start, end, idRegion = dictRegion[chrmName][i]
			startBP = random.randint(0, bpTotal)
			if start - startBP < chrms[chrmName][0]:
				startBP = chrms[chrmName][0] - start
				start = chrms[chrmName][0]
			else:
				start = start - startBP
			if end + (bpTotal - startBP) > chrms[chrmName][1]:
				startBP = end - chrms[chrmName][1]
				end = chrms[chrmName][1]
				if start - startBP >= chrms[chrmName][0]:
					start = start - startBP
			else:
				end = end + (bpTotal - startBP)
			F.write('%s\t%d\t%d\t%s\n' %(chrmName, start, end, idRegion))
			dictRegion[chrmName][i] = [start, end, idRegion]
	F.close()
