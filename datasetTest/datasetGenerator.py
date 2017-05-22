#!/usr/bin/python
# -*- coding: utf-8 -*-
import random, re, pprint


#################
### Functions ### ------------------------------------------------------------------------------------
#################

# Create random startCoord for regionsA without overlap bewteen regions
def randomStartA(regions, bpRegionA, bpRegionB, overlapSize, chrmStart, chrmEnd):
	# If no correct start position after 50 try, return True to change the chrm
	for nbTry in range(1, 10):
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
def accurateStartB(startA, bpRegionA, bpRegionB, overlapSize):
	if random.choice([True, False]):
		startCoord = startA - bpRegionB + overlapSize
	else:
		startCoord = startA + bpRegionA - overlapSize
	return startCoord

#
def randomStartB(startListA, startListB, bpRegionA, bpRegionB, overlapSize, chrmStart, chrmEnd):
	# If no correct start position after 50 try, return True to change the chrm
	for nbTry in range(1, 10):
		startCoord = random.randint(chrmStart, chrmEnd)
		notIntersect = True
		# Test if position don't overlap others regionA + possible regionB	
		for startRegion in startListA:	
			if startCoord >= startRegion-bpRegionB+overlapSize and startCoord <= startRegion+bpRegionA+bpRegionB-overlapSize:
				notIntersect = False
		for startRegion in startListB:	
			if startCoord >= startRegion and startCoord <= startRegion+bpRegionB:
				notIntersect = False
		if notIntersect:
			return startCoord
	return False


# Distribute a equal number of bp remaining across the regions, until bpRemaining=0
def getDistribution(bpRemaining, nbRegion):
	distribution = []
	for i in range(nbRegion):
		if bpRemaining != 0:
			bpRemaining -= 1
			distribution.append(1)
		else:
			distribution.append(0)
	while bpRemaining != 0:
		for i in range(nbRegion):
			if bpRemaining != 0:
				bpRemaining -= 1
				distribution[i] += 1
			else:
				break
	return distribution

#######################
### Input Variables ### ------------------------------------------------------------------------------------
#######################
genomePath = "./"
genomeName = "genome_test.bed"
# 1. Percentage of BP overlap between A and B
overlapPerc = 1.0
# 2. BP total for all regions (bpTotal >= nbRegion ; MAX = genomeSize/10000)
bpTotalA = 20000000
bpTotalB = 20000000
# 3. Number of regions 		Max = (bpTotalA/2)*overlapPerc < nbRegionA
nbRegionA = 10
nbRegionB = 100

folderPath = './'
nameA = 'datasetA_%d_%d.bed' %(nbRegionA, nbRegionB)
nameB = 'datasetB_%d_%d.bed' %(nbRegionA, nbRegionB)

# Quit if not enough bpTotalA or too much regions
if ((bpTotalA+bpTotalB)/2)*(overlapPerc/100) < 1:
	print "Error, not enough bpTotalA or too much regions!"
	exit()

############################
### Calculated Variables ### ------------------------------------------------------------------------------------
############################

# Calculate the number of BP by region
bpRegionA = bpTotalA // nbRegionA
bpRegionB = bpTotalB // nbRegionB

# Number of BP remaining after divided bpTotal by the number of region
bpRegionA_Remaining = bpTotalA - bpRegionA * nbRegionA
bpRegionB_Remaining = bpTotalB - bpRegionB * nbRegionB

# Distribute the bpRegion_Remaining to each regions
distriRegionA = getDistribution( bpRegionA_Remaining, nbRegionA )
distriRegionB = getDistribution( bpRegionB_Remaining, nbRegionB )

# Calculate the total number BP for all overlap
bpOverlapTotal = ((bpTotalA+bpTotalB)*(overlapPerc/100)) / 2
averageOverlap = bpOverlapTotal // min(nbRegionA,nbRegionB)

# Number of BP remaining after each regions have a least an averageOverlap bp
bpOverlap_Remaining = bpOverlapTotal - averageOverlap * min(nbRegionA,nbRegionB)

# Distribute the bpOverlap_Remaining to each regions
distriOverlap = getDistribution( bpOverlap_Remaining, min(nbRegionA,nbRegionB) )

# Calculate the maximum possible size of the regions and for the overlap
maxRegionA = bpRegionA + max(distriRegionA)
maxRegionB = bpRegionB + max(distriRegionB)
overlapSizeMax = averageOverlap + max(distriOverlap)

########################
### Genome File Part ### ------------------------------------------------------------------------------------
########################

# Save chroms coordinates and name of the genome file in a list
chrms = []
try:
	FG = open(genomePath+genomeName, 'r')
	for line in FG:
		if re.search('^track', line) is None:
			line = line.rstrip('\n')
			name, start, end = re.split('\t', line)
			if int(end)-int(start)+1 >= bpRegionA+bpRegionB*2+bpRegionA:
				chrms.append( [name, int(start), int(end)] )
except IOError:
    	print "Could not read file:", genomeName
    	exit()
finally:
   	FG.close()

#########################
### Main Part - FileA ### ------------------------------------------------------------------------------------
#########################

dictRegionA, dictRegionA_save = {}, {}
# Open FileA + write header
FA = open(folderPath+nameA, 'w')
FA.write('track	type=Bed\tname="%s"\n' % nameA)
# Create datasetA with regionsA
for idRegionA in range(nbRegionA):
	# Create region name and add a number of bp remaining for the region bp
	dataName = 'DataA_%d' %idRegionA
	completeRegionA = bpRegionA + distriRegionA[idRegionA]
	# Loop until find a start coordinate
	while True:
		# Get a random chromosom and this coordinates
		randomChrm = random.randint(0, len(chrms)-1)
		chrmName, chrmStart, chrmEnd = chrms[randomChrm]
		# Get the former start for this random chromosom
		startList = []
		if chrmName in dictRegionA:
			for e in dictRegionA[chrmName]:
				startList.append(e[0])
		# Generate a random start for this new region on the random chromosom
		startCoord = randomStartA(startList, completeRegionA, maxRegionB, overlapSizeMax, chrmStart, chrmEnd)
		if startCoord != False:
			break
	# Save this new region, the order of generating and the distribution
	if chrmName not in dictRegionA:
		dictRegionA[chrmName] = [ [startCoord, idRegionA, distriRegionA[idRegionA]] ]
		dictRegionA_save[chrmName] = [startCoord]
	else:
		dictRegionA[chrmName].append( [startCoord, idRegionA, distriRegionA[idRegionA]] )
		dictRegionA_save[chrmName].append( startCoord )
	# Calculate the coordinate of the end of the region
	endCoord = startCoord + completeRegionA - 1
	# Write the region into the A file
	FA.write('%s\t%d\t%d\t%s\n' %(chrmName, startCoord, endCoord, dataName))
FA.close()

#########################
### Main Part - FileB ### ------------------------------------------------------------------------------------
#########################

total = 0
dictRegionB = {}
minRegionAB = min(nbRegionA,nbRegionB)
remainingRegionB = nbRegionB - minRegionAB
# Open FileB + write header
FB = open(folderPath+nameB, 'w')
FB.write('track	type=Bed\tname="%s"\n' % nameB)
# Create datasetB with regionsB
for regionB in range(minRegionAB):
	while True:
		chrmName = random.choice(dictRegionA.keys())
		if len(dictRegionA[chrmName]) != 0:
			startCoordA, idRegionA, bpRegionADistri = dictRegionA[chrmName].pop(0)
			break
	dataName = 'DataB_%d' %idRegionA
	completeRegionA = bpRegionA + bpRegionADistri
	completeRegionB = bpRegionB + distriRegionB[regionB]
	overlapSize = averageOverlap + distriOverlap.pop(0)
	startCoordB = accurateStartB(startCoordA, completeRegionA, completeRegionB, overlapSize)
	endCoordB = startCoordB + completeRegionB - 1
	FB.write('%s\t%d\t%d\t%s\n' %(chrmName, startCoordB, endCoordB, dataName))
	if chrmName not in dictRegionB:
		dictRegionB[chrmName] = [startCoordB]
	else:
		dictRegionB[chrmName].append(startCoordB)

# Create regionB remaining with no overlap
for regionB in range(minRegionAB, minRegionAB+remainingRegionB):
	dataName = 'noOverlapB_%d' %regionB
	completeRegionB = bpRegionB + distriRegionB[regionB]
	# Loop until find a start coordinate
	while True:
		# Get a random chromosom and this coordinates
		randomChrm = random.randint(0, len(chrms)-1)
		chrmName, chrmStart, chrmEnd = chrms[randomChrm]
		# Get the former start for this random chromosom
		startListA, startListB = [], []
		if chrmName in dictRegionA_save:
			startListA = dictRegionA_save[chrmName]
		if chrmName in dictRegionB:
			startListB = dictRegionB[chrmName]
		# Generate a random start for this new region on the random chromosom
		startCoord = randomStartB(startListA, startListB, maxRegionA, maxRegionB, overlapSizeMax, chrmStart, chrmEnd)
		if startCoord != False:
			break
	# Save this new region, the order of generating and the distribution
	if chrmName not in dictRegionB:
		dictRegionB[chrmName] = [startCoord]
	else:
		dictRegionB[chrmName].append(startCoord)
	# Calculate the coordinate of the end of the region
	endCoord = startCoord + completeRegionB - 1
	FB.write('%s\t%d\t%d\t%s\n' %(chrmName, startCoord, endCoord, dataName))
FB.close()

