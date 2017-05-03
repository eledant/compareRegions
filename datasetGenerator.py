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
def randomStartB(regionA, genomeSize, bpRegion, bpOverlap):
	startCoord = regionA - bpRegion + bpOverlap
	if startCoord < 0:
		startCoord = regionA + bpRegion - bpOverlap
	return startCoord

#################
### Variables ###
#################
folderPath = './'
nameA = 'datasetA.txt'
nameB = 'datasetB.txt'

genomeSize = 100
# 1. BP overlap between A and B (condition = 0.1%, 1%, 10% ...)
bpOverlap = 4
# 2. BP for a region (<= bpOverlap, condition = same or different)
bpRegionA = 9
bpRegionB = 9
# 3. Number of regions (condition = same or different)
nbRegionA = 5
nbRegionB = 5


#################
### Main part ###
#################

# Create datasetA with regionsA
FH = open(folderPath+nameA, 'w')
FH.write('track	type=Bed\tname="%s"\n' % nameA)
startRegionA = []
for regionA in range(nbRegionA):
	regionName = 'Data%d' %regionA
	startCoord = randomStartA(startRegionA, genomeSize, bpRegionA)
	startRegionA.append(startCoord)
	endCoord = startCoord + bpRegionA - 1
	FH.write('chrm\t%d\t%d\t%s\n' %(startCoord, endCoord, regionName))
FH.close()

# Create datasetB with regionsB
FH = open(folderPath+nameB, 'w')
FH.write('track	type=Bed\tname="%s"\n' % nameB)
for regionB in range(nbRegionB):
	regionName = 'Data%d' %regionB
	startCoord = randomStartB(startRegionA.pop(0), genomeSize, bpRegionB, bpOverlap)
	endCoord = startCoord + bpRegionB - 1
	FH.write('chrm\t%d\t%d\t%s\n' %(startCoord, endCoord, regionName))
FH.close()



