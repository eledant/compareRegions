#!/usr/bin/python
# -*- coding: utf-8 -*-
import random, re, math

importPath = './'
importFileA = 'datasetA'
importFileB = 'datasetB'
genomePath = '../../../../'
genomeName = 'genome_structure'
chrms = {}
xTime = 10.0


# Save chroms coordinates and name of the genome file in a list
try:
	FG = open(genomePath+genomeName+'.bed', 'r')
	for line in FG:
		if re.search('^track', line) is None:
			line = line.rstrip('\n')
			name, start, end = re.split('\t', line)
			chrms[name] = [int(start), int(end)]
		else:
			firstLine = line
except IOError:
    	print "Could not read file:", genomeName
    	exit()
finally:
   	FG.close()


for x in range(0,int(xTime+1)):
	dictRegionA, dictRegionB = {}, {}
	# Save regions of an import fileA
	F = open(importPath+importFileA+'_'+str(x)+'.bed', 'r')
	for line in F:
		if re.search('^track', line) is None:
			line = line.rstrip('\n')
			chrmName, start, end, idRegion = re.split('\t', line)
			start, end, idRegion = int(start), int(end), idRegion
			if chrmName not in dictRegionA:
				dictRegionA[chrmName] = [ [start, end, idRegion] ]
			else:
				dictRegionA[chrmName].append( [start, end, idRegion] )
	F.close()

	# Save regions of an import fileB
	F = open(importPath+importFileB+'_'+str(x)+'.bed', 'r')
	for line in F:
		if re.search('^track', line) is None:
			line = line.rstrip('\n')
			chrmName, start, end, idRegion = re.split('\t', line)
			start, end, idRegion = int(start), int(end), idRegion
			if chrmName not in dictRegionB:
				dictRegionB[chrmName] = [ [start, end, idRegion] ]
			else:
				dictRegionB[chrmName].append( [start, end, idRegion] )
	F.close()

	# Sort regions
	for chrmName in dictRegionA:
		dictRegionA[chrmName].sort(key=lambda x: x[0])
		dictRegionB[chrmName].sort(key=lambda x: x[0])

	# Create new fileA and fileB
	F = open(genomeName+'_'+str(x)+'.bed', 'w')
	F.write(firstLine)
	for chrmName in dictRegionA:
		startChrm = min(dictRegionA[chrmName][0][0], dictRegionB[chrmName][0][0])
		endChrm = max(dictRegionA[chrmName][-1][1], dictRegionB[chrmName][-1][1])
		F.write('%s\t%d\t%d\n' %(chrmName+'_begin', chrms[chrmName][0], startChrm-1))
		F.write('%s\t%d\t%d\n' %(chrmName, startChrm, endChrm))
		F.write('%s\t%d\t%d\n' %(chrmName+'_end', endChrm+1, chrms[chrmName][1]))
	F.close()



