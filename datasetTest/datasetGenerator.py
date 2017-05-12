#!/usr/bin/python
# -*- coding: utf-8 -*-
import random, re, pprint


#################
### Functions ###
#################

# Create random startCoord for regionsA without overlap bewteen regions
def randomStartA(regions, bpRegionA, bpRegionB, overlapSize, chrmStart, chrmEnd):
	# If no correct start position after 50 try, return True to change the chrm
	for nbTry in range(1, 50):
		startCoord = random.randint( chrmStart+bpRegionB, chrmEnd-bpRegionB )
		if not regions:
			return startCoord
		notIntersect = True
		# Test if position don't overlap others regionA + possible regionB	
		for startRegion in regions:	
			if startCoord+bpRegionA+bpRegionB-overlapSize >= startRegion-bpRegionB+overlapSize and startCoord-bpRegionB+overlapSize <= startRegion+bpRegionA+bpRegionB-overlapSize:
				notIntersect = False
		if notIntersect:
			return startCoord
	return False

# Create accurate startCoord for regionB, overlaping with regionA
def accurateStartB(regionA, bpRegion, overlapSize):
	if random.choice([True, False]):
		startCoord = regionA - bpRegion + overlapSize
	else:
		startCoord = regionA + bpRegion - overlapSize
	return startCoord

#################
### Variables ###
#################
genomePath = "./"
genomeName = "genome_test.bed"
# 1. Percentage of BP overlap between A and B
overlapPerc = 1.0
# 2. BP total for all regions (bpTotal >= nbRegion ; MAX = genomeSize/10000)
bpTotalA = 200000
bpTotalB = 200000
# 3. Number of regions 		Max = (bpTotalA/2)*overlapPerc < nbRegionA
nbRegionA = 50
nbRegionB = 50

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

# Save chrom coordinates of the genome
chrms = []
try:
	FG = open(genomePath+genomeName, 'r')
	for line in FG:
		if re.search('^track', line) is None:
			line = line.rstrip('\n')
			name, start, end = re.split('\t', line)
			if int(end)-int(start)+1 >= bpRegionA+bpRegionB*2-overlapSize:
				chrms.append( [name, int(start), int(end)] )
except IOError:
    	print "Could not read file:", genomeName
    	exit()
finally:
   	FG.close()

# Create datasetA with regionsA
regionStarts, regionOrders = {}, {}
FA = open(folderPath+nameA, 'w')
FA.write('track	type=Bed\tname="%s"\n' % nameA)
for regionA in range(nbRegionA):
	dataName = 'DataA_%d' %regionA
	# Loop until find a start coordinate
	while True:
		# Get a random chromosom and this coordinates
		randomChrm = random.randint(0, len(chrms)-1)
		chrmName, chrmStart, chrmEnd = chrms[randomChrm]
		# Calculate a random start for a new region on the chromosom
		if chrmName in regionStarts:
			regionStart = regionStarts[chrmName]
		else:
			regionStart = []
		startCoord = randomStartA(regionStart, bpRegionA, bpRegionB, overlapSize, chrmStart, chrmEnd)
		if startCoord != False:
			break
	# Save this new region
	if chrmName not in regionStarts:
		regionStarts[chrmName] = [startCoord]
		regionOrders[chrmName] = [regionA]
	else:
		regionStarts[chrmName].append(startCoord)
		regionOrders[chrmName].append(regionA)
	# Calculate the coordinate of the end of the region
	endCoord = startCoord + bpRegionA - 1
	# Write the region into the A file
	FA.write('%s\t%d\t%d\t%s\n' %(chrmName, startCoord, endCoord, dataName))
FA.close()

# Create datasetB with regionsB
regionB = 0
FB = open(folderPath+nameB, 'w')
FB.write('track	type=Bed\tname="%s"\n' % nameB)
for chrmName in regionStarts:
	for startA in regionStarts[chrmName]:
		dataName = 'DataB_%d' %regionOrders[chrmName].pop(0)
		startCoord = accurateStartB(startA, bpRegionB, overlapSize)
		endCoord = startCoord + bpRegionB - 1
		FB.write('%s\t%d\t%d\t%s\n' %(chrmName, startCoord, endCoord, dataName))
		regionB += 1
FB.close()



