#!/usr/bin/python
# -*- coding: utf-8 -*-
import random, re, math

importPath = './'
importFileA = 'datasetA'
importFileB = 'datasetB'
genomePath = '../../../../'
genomeName = 'genome_structure'
chrms, score = {}, {}
startA, endA = {}, {}
xTime = 10.0


# Save chroms coordinates and name of the genome file in a list
try:
	FG = open(genomePath+genomeName+'.bed', 'r')
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
		else:
			firstLineA = line
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
		else:
			firstLineB = line
	F.close()

	# Sort regions
	for chrmName in dictRegionA:
		dictRegionA[chrmName].sort(key=lambda x: x[0])
		dictRegionB[chrmName].sort(key=lambda x: x[0])

	# Calulate mean score
	for chrmName in dictRegionA:
		score[chrmName] = 0
		for region in dictRegionA[chrmName]:
			score[chrmName] += region[1]-region[0]+1

	# Create new fileA and fileB
	FA = open(importPath+'datasetOneA_'+str(x)+'.bed', 'w')
	FB = open(importPath+'datasetOneB_'+str(x)+'.bed', 'w')
	FA.write(firstLineA)
	FB.write(firstLineB)
	for chrmName in dictRegionA:
		start, end = min(dictRegionA[chrmName][0][0], dictRegionB[chrmName][0][0]), max(dictRegionA[chrmName][-1][1], dictRegionB[chrmName][-1][1])
		idRegionA, idRegionB = 'regionA_'+chrmName+'_'+str(x), 'regionB_'+chrmName+'_'+str(x)
		size = (end - start + 1)
		scoreBP = math.sqrt( score[chrmName] / float(size) )
		# Print overlaping regions
		FA.write('%s\t%d\t%d\t%s\t%g\n' %(chrmName, start, end, idRegionA, scoreBP ))
		FB.write('%s\t%d\t%d\t%s\t%g\n' %(chrmName, start, end, idRegionB, scoreBP ))
	FA.close()
	FB.close()

