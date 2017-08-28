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
		startCoord = random.randint( chrmStart+bpRegionB, chrmEnd-bpRegionA-bpRegionB )
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
	#if bpRegionA-bpRegionB > 0:
	#	averageMiddle = random.randint( 0, bpRegionA-bpRegionB )
	#else:
	averageMiddle = 0
	if random.choice([True, False]):
		startCoord = startA - bpRegionB + overlapSize + averageMiddle
	else:
		startCoord = startA + bpRegionA - overlapSize - averageMiddle
	return startCoord

# Create a new regionB which don't overlpa other regions
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
def getDistributionEqual(bpRemaining, nbRegion):
	distribution = []
	for i in range(nbRegion):
		if bpRemaining != 0:
			bpRemaining -= 1
			distribution.append(1)
		else:
			distribution.append(0)
	return distribution

# Distribute a random number of bp remaining across the regions, until bpRemaining=0
def getDistributionRandom1(bpRemaining, nbRegion, averageOverlap):
	distribution = []
	for i in range(nbRegion):
		if bpRemaining != 0:
			if bpRemaining < averageOverlap:
				bp = random.randint(0, bpRemaining)
			else:
				bp = random.randint(0, averageOverlap)
			bpRemaining -= bp
			distribution.append(bp)
		else:
			distribution.append(0)
	while bpRemaining != 0:
		for i in range(nbRegion):
			if bpRemaining == 0:
				break
			if bpRemaining < averageOverlap:
				bp = random.randint(0, bpRemaining)
			else:
				bp = random.randint(0, averageOverlap)
			bpRemaining -= bp
			distribution[i] += bp
	return distribution

# Distribute a random number of bp based on the regions size, until bpRemaining=0
def getDistributionRandom2(bpRemaining, nbRegion, bpRegion, distriRegion):
	distribution = []
	for i in range(nbRegion):
		regionSize = bpRegion + distriRegion[i]
		if bpRemaining != 0:
			if bpRemaining < regionSize:
				bp = random.randint(0, bpRemaining)
			else:
				bp = random.randint(0, regionSize)
			bpRemaining -= bp
			distribution.append(bp)
		else:
			distribution.append(0)
	while bpRemaining != 0:
		for i in range(nbRegion):
			regionSize = bpRegion + distriRegion[i]
			if bpRemaining == 0:
				break
			if bpRemaining < regionSize:
				bp = random.randint(0, bpRemaining)
			else:
				bp = random.randint(0, regionSize)
			bpRemaining -= bp
			distribution[i] += bp
	return distribution

# Distribute a random number of bp based on the number of regions, until bpRemaining=0
def getDistributionRandom3(bpRemaining, nbRegion):
	distribution = []
	for i in range(nbRegion):
		distribution.append(1)
	bpRemaining -= nbRegion
	while bpRemaining != 0:
		for i in range(nbRegion):
			if bpRemaining == 0:
				break
			bp = random.randint(0, bpRemaining)
			bpRemaining -= bp
			distribution[i] += bp
	return distribution

#######################
### Input Variables ### ------------------------------------------------------------------------------------
#######################
genomePath = "./"
genomeName = "genome_structure.bed"
# 1. Percentage of BP overlap between A and B
overlapPerc = 100.0
# 2. BP total for all regions (bpTotal >= nbRegion ; MAX = genomeSize/10000)
bpTotalA = 500000
bpTotalB = 500000
# 3. Number of regions 		Max = (bpTotalA/2)*overlapPerc < nbRegionA	
nbRegionA = 5000
nbRegionB = 5000

# 4. Overlap Distribution Type | 1=Equal, 2=AverageRandom, 3=TotalRandom
distriType = 1

# 5. BP Total Overlap Calculation | 1=Average, 2=sizeA, 3=sizeB
overlapCalcul = 1

# 6. BP by region Distribution | 1=Equal, 2=AverageRandom, 3=TotalRandom
bpRegionCalcul = 1

# Work with a existing datasetA		#Test4 -> start with low values
importA = False
importPath = './'
importFile = 'datasetA.bed'

folderPath = './'
nameA = 'datasetA_%s.bed' %(str(overlapPerc))
if not importA:
	nameB = 'datasetB_%s.bed' %(str(overlapPerc))
else:
	nameB = 'datasetB_%d.bed' %(overlapPerc)


# Quit if not enough bpTotalA or too much regions
if ((bpTotalA+bpTotalB)/2)*(overlapPerc/100) < 1:
	print "Error, not enough bpTotalA or too much regions!"

############################
### Calculated Variables ### ------------------------------------------------------------------------------------
############################

# Equal distribution of BP for each region
if bpRegionCalcul == 1:
	# Calculate the number of BP by region
	bpRegionA = bpTotalA // nbRegionA
	bpRegionB = bpTotalB // nbRegionB

	# Number of BP remaining after divided bpTotal by the number of region
	bpRegionA_Remaining = bpTotalA - bpRegionA * nbRegionA
	bpRegionB_Remaining = bpTotalB - bpRegionB * nbRegionB

	# Distribute the bpRegion_Remaining to each regions
	distriRegionA = getDistributionEqual( bpRegionA_Remaining, nbRegionA )
	distriRegionB = getDistributionEqual( bpRegionB_Remaining, nbRegionB )
# Total Random distribution of BP for each region
elif bpRegionCalcul == 2:
	distriRegionA = getDistributionRandom1( bpTotalA, nbRegionA, bpTotalA//nbRegionA )
	distriRegionB = getDistributionRandom1( bpTotalB, nbRegionB, bpTotalB//nbRegionB )
	bpRegionA = 0
	bpRegionB = 0

# Total Random distribution of BP for each region
elif bpRegionCalcul == 3:
	distriRegionA = getDistributionRandom3( bpTotalA, nbRegionA )
	distriRegionB = getDistributionRandom3( bpTotalB, nbRegionB )
	bpRegionA = 0
	bpRegionB = 0

# Calculate the total number BP for all overlap according to different methods
if overlapCalcul == 1:
	bpOverlapTotal = ((bpTotalA+bpTotalB)/2)*(overlapPerc/100)
elif overlapCalcul == 2:
	bpOverlapTotal = bpTotalA*(overlapPerc/100)
elif overlapCalcul == 3:
	bpOverlapTotal = bpTotalB*(overlapPerc/100)
# Calculate an average overlap for each region
averageOverlap = bpOverlapTotal // min(nbRegionA,nbRegionB)

# Equal number of overlap bp
if distriType == 1:
	# Number of BP remaining after each regions have a least an averageOverlap bp
	bpOverlap_Remaining = bpOverlapTotal - averageOverlap * min(nbRegionA,nbRegionB)

	# Distribute the bpOverlap_Remaining to each regions
	distriOverlap = getDistributionEqual( bpOverlap_Remaining, min(nbRegionA,nbRegionB) )
# Random number of bp based on the average overlap
elif distriType == 2:
	distriOverlap = getDistributionRandom1( bpOverlapTotal, min(nbRegionA,nbRegionB), averageOverlap )
	averageOverlap = 0

# Random number of bp based on the region size
elif distriType == 3:
	if bpRegionA+max(distriRegionA) < bpRegionB+max(distriRegionB):
		distriOverlap = getDistributionRandom2( bpOverlapTotal, min(nbRegionA,nbRegionB), bpRegionA,  distriRegionA )
	else:
		distriOverlap = getDistributionRandom2( bpOverlapTotal, min(nbRegionA,nbRegionB), bpRegionB,  distriRegionB )
	averageOverlap = 0

# Calculate the maximum possible size of the regions and for the overlap
maxRegionA = bpRegionA + max(distriRegionA)
maxRegionB = bpRegionB + max(distriRegionB)
overlapSizeMin = averageOverlap + min(distriOverlap)


########################
### Genome File Part ### ------------------------------------------------------------------------------------
########################

# Save chroms coordinates and name of the genome file in a list
chrms = []
coordRemap = 0
try:
	FG = open(genomePath+genomeName, 'r')
	for line in FG:
		if re.search('^track', line) is None:
			line = line.rstrip('\n')
			name, start, end = re.split('\t', line)
			if int(end)-int(start)+1 >= ( (bpRegionA+max(distriRegionA))+(bpRegionB+max(distriRegionB)) )*2:		
				chrms.append( [name, int(start), int(end), coordRemap] )
				coordRemap += int(end) - int(start)
except IOError:
    	print "Could not read file:", genomeName
    	exit()
finally:
   	FG.close()



#########################
### Main Part - FileA ### ------------------------------------------------------------------------------------
#########################

dictRegionA, dictRegionA_save = {}, {}

# Save regions of an import fileA
if importA:
	# The number of bp by region is now in 'dictRegionA'
	bpRegionA = 0
	nbRegionA = 0
	FA = open(importPath+importFile, 'r')
	for line in FA:
		if re.search('^track', line) is None:
			nbRegionA += 1
			line = line.rstrip('\n')
			chrmName, start, end, idRegion = re.split('\t', line)
			idRegion = re.split('_', idRegion)[-1]
			start, end, idRegion = int(start), int(end), int(idRegion)
			if chrmName not in dictRegionA:
				dictRegionA[chrmName] = [ [start, idRegion, end-start+1] ]
				dictRegionA_save[chrmName] = [start]
			else:
				dictRegionA[chrmName].append( [start, idRegion, end-start+1] )
				dictRegionA_save[chrmName].append( start )
	FA.close()
	


# Generate a fileA from variables
else:
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
			randomStart = random.randint(0, coordRemap)
			for chrmInfo in chrms:
				chrmName, chrmStart, chrmEnd, startCoord = chrmInfo
				if randomStart > startCoord and randomStart < (startCoord+chrmEnd-chrmStart):
					break
			# Get the former start for this random chromosom
			startList = []
			if chrmName in dictRegionA:
				for e in dictRegionA[chrmName]:
					startList.append(e[0])
			# Generate a random start for this new region on the random chromosom
			startCoord = randomStartA(startList, completeRegionA, maxRegionB, overlapSizeMin, chrmStart, chrmEnd)
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

dictRegionB = {}
minRegionAB = min(nbRegionA,nbRegionB)
remainingRegionB = nbRegionB - minRegionAB
# Open FileB + write header
FB = open(folderPath+nameB, 'w')
FB.write('track	type=Bed\tname="%s"\n' % nameB)
# Create datasetB with regionsB
for regionB in range(minRegionAB):
	# Take a random A region
	if bpRegionCalcul == 1:
		while True:
			chrmName = random.choice(dictRegionA.keys())
			if len(dictRegionA[chrmName]) != 0:
				startCoordA, idRegionA, bpRegionADistri = dictRegionA[chrmName].pop(0)
				break
	# Special case for bpRegionCalcul method 2 and 3
	elif bpRegionCalcul == 2 or bpRegionCalcul == 3:
		maxBP = 0
		for key in dictRegionA:
			for i in range(len(dictRegionA[key])):
				if dictRegionA[key][i][2] > maxBP:
					maxBP = dictRegionA[key][i][2]
					chrmName = key
					position = i
		startCoordA, idRegionA, bpRegionADistri = dictRegionA[chrmName].pop(position)		 
	# Calculate correct size of A region
	completeRegionA = bpRegionA + bpRegionADistri
	# Calculate correct size of B region
	if bpRegionCalcul == 1:
		completeRegionB = bpRegionB + distriRegionB[regionB]
	elif bpRegionCalcul == 2 or bpRegionCalcul == 3:
		completeRegionB = bpRegionB + max(distriRegionB)
		distriRegionB.remove(max(distriRegionB))
	dataName = 'DataB_%d' %idRegionA
	# Take a random overlap
	if bpRegionCalcul == 1:
		overlapSize = averageOverlap + distriOverlap.pop(0)
	elif bpRegionCalcul == 2 or bpRegionCalcul == 3:
		overlapSize = averageOverlap + max(distriOverlap)
		distriOverlap.remove(max(distriOverlap))
	# Calculate coordinates
	startCoordB = accurateStartB(startCoordA, completeRegionA, completeRegionB, overlapSize)
	endCoordB = startCoordB + completeRegionB - 1
	FB.write('%s\t%d\t%d\t%s\n' %(chrmName, startCoordB, endCoordB, dataName))
	if chrmName not in dictRegionB:
		dictRegionB[chrmName] = [startCoordB]
	else:
		dictRegionB[chrmName].append(startCoordB)

exit()

# Create regionB remaining with no overlap
for regionB in range(minRegionAB, minRegionAB+remainingRegionB):
	dataName = 'noOverlapB_%d' %regionB
	if bpRegionCalcul == 1:
		completeRegionB = bpRegionB + distriRegionB[regionB]
	elif bpRegionCalcul == 2 or bpRegionCalcul == 3:
		completeRegionB = bpRegionB + max(distriRegionB)
		distriRegionB.remove(max(distriRegionB))
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
		startCoord = randomStartB(startListA, startListB, maxRegionA, maxRegionB, overlapSizeMin, chrmStart, chrmEnd)
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

