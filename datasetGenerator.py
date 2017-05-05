#!/usr/bin/python
# -*- coding: utf-8 -*-
import random


#################
### Functions ###
#################

# Create random startCoord for regionsA without overlap bewteen regions
def randomStartA(regions, genomeSize, bpRegion):
	while True:
		startCoord = random.randint( 0, genomeSize-bpRegion )
		if not regions:
			return startCoord
		for startRegion in regions:	
			if startCoord < startRegion or startCoord > startRegion+bpRegion-1:
				return startCoord

# Create accurate startCoord for regionB, overlaping with regionA
def accurateStartB(regionA, genomeSize, bpRegion, overlapSize):
	startCoord = regionA - bpRegion + overlapSize
	if startCoord < 0:
		startCoord = regionA + bpRegion - overlapSize
	return startCoord

#################
### Variables ###
#################
genomeSize = 2600000000
# 1. Percentage of BP overlap between A and B
overlapPerc = 10.0
# 2. BP total for all regions (bpTotal >= nbRegion)
bpTotalA = 100000000
bpTotalB = 100000000
# 3. Number of regions 		Max = (bpTotalA/2)*overlapPerc < nbRegionA
nbRegionA = 100000
nbRegionB = 100000

folderPath = './'
nameA = 'datasetA_%d.bed' %nbRegionA
nameB = 'datasetB_%d.bed' %nbRegionB

#################
### Main part ###
#################

# Calculate the number of BP by region
bpRegionA = bpTotalA / nbRegionA
bpRegionB = bpTotalB / nbRegionB
# Calculate the BP of an overlap
bpOverlap = ((bpTotalA+bpTotalB)/2)*(overlapPerc/100)
overlapSize = bpOverlap / ((nbRegionA+nbRegionB)/2)

# Quit if not enough bpTotalA or too much regions
if (bpTotalA/2)*overlapPerc < nbRegionA or (bpTotalB/2)*overlapPerc < nbRegionB:
	print "Error, not enough bpTotalA or too much regions!"
	exit()

# Create datasetA with regionsA
FA = open(folderPath+nameA, 'w')
FA.write('track	type=Bed\tname="%s"\n' % nameA)
startRegionA = []
for regionA in range(nbRegionA):
	regionName = 'Data%d' %regionA
	startCoord = randomStartA(startRegionA, genomeSize, bpRegionA)
	startRegionA.append(startCoord)
	endCoord = startCoord + bpRegionA - 1
	FA.write('chrm\t%d\t%d\t%s\n' %(startCoord, endCoord, regionName))
FA.close()

# Create datasetB with regionsB
FB = open(folderPath+nameB, 'w')
FB.write('track	type=Bed\tname="%s"\n' % nameB)
for regionB in range(nbRegionB):
	regionName = 'Data%d' %regionB
	startCoord = accurateStartB(startRegionA.pop(0), genomeSize, bpRegionB, overlapSize)
	endCoord = startCoord + bpRegionB - 1
	FB.write('chrm\t%d\t%d\t%s\n' %(startCoord, endCoord, regionName))
FB.close()



