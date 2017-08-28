#!/usr/bin/python
# -*- coding: utf-8 -*-
import random, re, math

option = '_0'
importPath = './'
importFileA = 'datasetA'
importFileB = 'datasetB'
dictRegionA, dictRegionB, bary, distance = {}, {}, {}, {}
xTime = 10.0


# Save regions of an import fileA
F = open(importPath+importFileA+option+'.bed', 'r')
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
F = open(importPath+importFileB+option+'.bed', 'r')
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


# Sort regions and calculate barycenters
for chrmName in dictRegionA:
	dictRegionA[chrmName].sort(key=lambda x: x[0])
	dictRegionB[chrmName].sort(key=lambda x: x[0])
	bary[chrmName] = 0
	for i in range(len(dictRegionA[chrmName])):
		startA, endA, idA = dictRegionA[chrmName][i]
		startB, endB, idB = dictRegionB[chrmName][i]
		bary[chrmName] += min(startA, startB) + ( max(endA, endB)-min(startA, startB) ) / 2
	bary[chrmName] /= len(dictRegionA[chrmName])


partOne, partTwo = {}, {}
# Calculate distance barycenter/region and sort regions
for chrmName in dictRegionA:
	bar = bary[chrmName]
	distance[chrmName] = {}
	firstTime = True
	for i in range(len(dictRegionA[chrmName])):
		startA, endA, idA = dictRegionA[chrmName][i]
		startB, endB, idB = dictRegionB[chrmName][i]
		midPoint = min(startA, startB) + ( max(endA, endB)-min(startA, startB) ) / 2
		distance[chrmName][i] = int(bar - midPoint)
		if midPoint > bar and firstTime:
			firstTime = False
			if i == 0:
				partOne[chrmName] = 0
			else:
				partOne[chrmName] = i-1
			partTwo[chrmName] = i
	if firstTime:
		partOne[chrmName] = 0
		partTwo[chrmName] = 0


# Create new fileA and fileB
for x in range(1,int(xTime+1)):
	FA = open('datasetA_'+str(x)+'.bed', 'w')
	FB = open('datasetB_'+str(x)+'.bed', 'w')
	FA.write(firstLineA)
	FB.write(firstLineB)
	for chrmName in dictRegionA:
		prevOneStart, prevOneEnd = 0, 0
		prevTwoStart, prevTwoEnd = 0, 0
		firstOne = True
		for i in reversed(range(0, partOne[chrmName]+1)):
			dist = abs( int( (distance[chrmName][i] * x) / xTime ) )
			startA, endA, idA = dictRegionA[chrmName][i]
			startB, endB, idB = dictRegionB[chrmName][i]
			if max(endA, endB)+dist >= prevOneStart and not firstOne:
				if startA > startB:
					FA.write('%s\t%d\t%d\t%s\n' %(chrmName, prevOneStart-1-(endA-startA), prevOneStart-1, idA))
					FB.write('%s\t%d\t%d\t%s\n' %(chrmName, startB+(prevOneStart-1-endA), endB+(prevOneStart-1-endA), idB))
					prevOneStart, prevOneEnd = startB+(prevOneStart-1-endA), prevOneStart-1
				else:
					FA.write('%s\t%d\t%d\t%s\n' %(chrmName, startA+(prevOneStart-1-endB), endA+(prevOneStart-1-endB), idA))
					FB.write('%s\t%d\t%d\t%s\n' %(chrmName, prevOneStart-1-(endB-startB), prevOneStart-1, idB))
					prevOneStart, prevOneEnd = startA+(prevOneStart-1-endB), prevOneStart-1
			else:
				FA.write('%s\t%d\t%d\t%s\n' %(chrmName, startA+dist, endA+dist, idA))
				FB.write('%s\t%d\t%d\t%s\n' %(chrmName, startB+dist, endB+dist, idB))
				prevOneStart, prevOneEnd = min(startA, startB)+dist, max(endA, endB)+dist
				if firstOne:
					firstOne = False
					prevTwoStart, prevTwoEnd = prevOneStart, prevOneEnd
		for i in range(partTwo[chrmName], len(dictRegionA[chrmName])):
			dist = abs( int( (distance[chrmName][i] * x) / xTime ) )
			startA, endA, idA = dictRegionA[chrmName][i]
			startB, endB, idB = dictRegionB[chrmName][i]
			if min(startA, startB)-dist <= prevTwoEnd:
				if startA > startB:
					FA.write('%s\t%d\t%d\t%s\n' %(chrmName, startA-(startB-prevTwoEnd-1), endA-(startB-prevTwoEnd-1), idA))
					FB.write('%s\t%d\t%d\t%s\n' %(chrmName, prevTwoEnd+1, prevTwoEnd+1+(endB-startB), idB))
					prevTwoStart, prevTwoEnd = prevTwoEnd+1, endA-(startB-prevTwoEnd-1)
				else:
					FA.write('%s\t%d\t%d\t%s\n' %(chrmName, prevTwoEnd+1, prevTwoEnd+1+(endA-startA), idA))
					FB.write('%s\t%d\t%d\t%s\n' %(chrmName, startB-(startA-prevTwoEnd-1), endB-(startA-prevTwoEnd-1), idB))
					prevTwoStart, prevTwoEnd = prevTwoEnd+1, endB-(startA-prevTwoEnd-1)
			else:
				FA.write('%s\t%d\t%d\t%s\n' %(chrmName, startA-dist, endA-dist, idA))
				FB.write('%s\t%d\t%d\t%s\n' %(chrmName, startB-dist, endB-dist, idB))
				prevTwoStart, prevTwoEnd = min(startA, startB)-dist, max(endA, endB)-dist
	FA.close()
	FB.close()



