import rpy2.robjects as ro
import math
import numpy as np

# A list which contains the overlaps size and the split regions name
class OverlapMatrix(list):

	#############################
	### Initialize the matrix ###
	#############################
	def __init__(self, nA, nB, sizeA, sizeB):
		self.indexA_list, self.indexB_list = [], []
		self.splitRegionA, self.splitRegionB = [], []
		# Number of regions in datasetA and datasetB
		self.dataA_nb, self.dataB_nb = float(nA), float(nB)
		# Total size of all the regions for each dataset 
		self.dataA_size, self.dataB_size = float(sizeA), float(sizeB)
		# For Pearson calcul
		self.scoreA, self.scoreB = [], []


	###################################################
	### Add the size of a new overlap to the matrix ###
	###################################################
	def addOverlap(self, overlap, indexA, indexB):
		self.append(overlap)
		self.indexA_list.append(indexA)
		self.indexB_list.append(indexB)


	###################################
	### Add a list of split regions ###
	###################################
	def addSplitRegions(self, _splitRegionA, _splitRegionB):
		self.splitRegionA = list(set(_splitRegionA))
		self.splitRegionB = list(set(_splitRegionB))


	#############################
	### Set the overlaps size ###
	#############################
	def setOverlapMeasures(self):
		self.overlapSize, self.overlapsRegionNumber = 0, 0
		for i in range(len(self)):
			self.overlapSize += self[i]
			self.overlapsRegionNumber += 1
	

	######################################################
	### Return the number of A regions with an overlap ###
	######################################################
	def getOverA(self):
		return len(list(set(self.indexA_list))) - len(list(set(self.splitRegionA)))


	######################################################
	### Return the number of B regions with an overlap ###
	######################################################
	def getOverB(self):
		return len(list(set(self.indexB_list))) - len(list(set(self.splitRegionB)))
	

	############################################
	### Set the size of the reference genome ###
	############################################
	def addGenomeSize(self, _genomeSize):
		self.genomeSize = float(_genomeSize)

	
	########################################################
	### Get variables for the Pearson coefficient method ###
	########################################################
	def getPearsonVariables(self, regionsA, regionsB):
		self.pA, self.pB, self.pA_sq, self.pB_sq = 0, 0, 0, 0	
		for dataA in regionsA:
			scoreA = 1
			if 'score' in dataA:
				scoreA = dataA['score']
			self.pA += (dataA['endCoord']-dataA['startCoord']+1) * scoreA
			self.pA_sq += (dataA['endCoord']-dataA['startCoord']+1) * pow(scoreA,2)
		for dataB in regionsB:
			scoreB = 1
			if 'score' in dataB:
				scoreB = dataB['score']
			self.pB += (dataB['endCoord']-dataB['startCoord']+1) * scoreB
			self.pB_sq += (dataB['endCoord']-dataB['startCoord']+1) * pow(scoreB,2)

# -----------------------------------------------------------------------------------------

	########################################################################################
	### Menu : select the method to calculate the score according to the selected option ###
	########################################################################################
	def calcStats(self, randMatrices, fileB_name, args):
		fileB_name = fileB_name.split('/')[-1]
		# 0=line, 1=score_mean, 2=distribution, 3=expected_score, 4=title
		res_stats = []
		if args['-l'] in ['all', 'def']:
			res_stats.append( self.scoreProduct(randMatrices, fileB_name, args) )
		if args['-l'] in ['all', 'psn']:
			res_stats.append( self.pearson(randMatrices, fileB_name, args) )
		if args['-l'] in ['all', 'npm']:
			res_stats.append( self.npmi(randMatrices, fileB_name, args) )
		if args['-l'] in ['all', 'jac']:
			res_stats.append( self.jaccard(randMatrices, fileB_name, args) )
		if args['-l'] in ['all', 'ebo']:
			res_stats.append( self.encodeBaseOver(randMatrices, fileB_name, args) )
		if args['-l'] in ['all', 'ero']:
			res_stats.append( self.encodeRegionOver(randMatrices, fileB_name, args) )
		if args['-l'] in ['all', 'pwe']:
			res_stats.append( self.ernstKellis(randMatrices, fileB_name, args) )
		if args['-l'] in ['all', 'pbn']:
			res_stats.append( self.pBinom(randMatrices, fileB_name, args) )
		return res_stats

	################################################
	### Return the z-score, p-value, mean and sd ###
	################################################
	def getZScore(self, value, distribution):
		mean = np.mean(np.array(distribution))
		sd = np.std(np.array(distribution))
		if sd != 0:
			zscore = (value-mean)/sd
		else:
			zscore = 0
		return [zscore, mean, sd]

	
	######################################################
	### Calculate the z-score and return a output_line ###
	######################################################
	def scoreProduct(self, randMatrices, fileB_name, args):
		# Get the overlaps size and the number of overlap in A and B datasets
		randOverlap = []
		for m in randMatrices:
			randOverlap.append( m.overlapSize )
		zscore, mean, sd = self.getZScore(self.overlapSize, randOverlap)
		
		# Create output for the <B_file>
		if args['-l'] == 'all':
			output = "\tScore Product\t%g\t%g\t%g" %(self.overlapSize, zscore, mean)
		else:
			output = "#2: %.2f\t%.2f\t%d/%d\t%d: %d/%d\t%d: %d/%d\t%s" %(zscore, mean, self.overlapSize, self.genomeSize, self.dataA_size, self.getOverA(), self.dataA_nb, self.dataB_size,self.getOverB(),  self.dataB_nb, fileB_name)
		return [output, mean, randOverlap, self.overlapSize, 'Score Product', zscore, 0, min(self.dataA_size, self.dataB_size)]


	#################################################################
	### Calculate the Jaccard similarity and return a output_line ###
	#################################################################
	def jaccard(self, randMatrices, fileB_name, args):
		def calcJaccard(pAB, pA, pB):
			if pA + pB - pAB == 0:
				return 0
			else:
				return pAB / (pA + pB - pAB)
		
		jaccard = calcJaccard(self.overlapSize, self.dataA_size, self.dataB_size)
		jaccard_l = []
		for m in randMatrices:
			jaccard_l.append( calcJaccard(m.overlapSize, m.dataA_size, m.dataB_size) )
		jaccard_mean = sum(jaccard_l) / len(jaccard_l)
		zscore = self.getZScore(jaccard, jaccard_l)[0]
		jaccard_min = calcJaccard(0, self.dataA_size, self.dataB_size)
		jaccard_max = calcJaccard(min(self.dataA_size, self.dataB_size), self.dataA_size, self.dataB_size)

		if args['-l'] == 'all':
			output = "\tJaccard\t%g\t%g\t%g" %(jaccard, zscore, jaccard_mean)
		else:
			output = "#2: %.2f\t%.2f\t%d/%d\t%d: %d/%d\t%d: %d/%d\t%s" %(zscore, jaccard_mean, jaccard, self.genomeSize, self.dataA_size, self.getOverA(), self.dataA_nb, self.dataB_size,self.getOverB(),  self.dataB_nb, fileB_name)
		return [output, jaccard_mean, jaccard_l, jaccard, 'Jaccard', zscore, jaccard_min, jaccard_max]


	##########################################################################
	### Calculate the basepair overlap statistics and return a output_line ###
	##########################################################################
	def encodeBaseOver(self, randMatrices, fileB_name, args):
		def calcBaseOver(pAB, pA):
			return pAB / pA

		baseOver = calcBaseOver(self.overlapSize, self.dataA_size)
		baseOver_l = []
		for m in randMatrices:
			baseOver_l.append( calcBaseOver( m.overlapSize, m.dataA_size) )
		baseOver_mean = sum(baseOver_l) / len(baseOver_l)
		zscore = self.getZScore(baseOver, baseOver_l)[0]
		baseOver_min = calcBaseOver(0, self.dataA_size)
		baseOver_max = calcBaseOver(min(self.dataA_size, self.dataB_size), self.dataA_size)

		if args['-l'] == 'all':
			output = "\tBases\t%g\t%g\t%g" %(baseOver, zscore, baseOver_mean)
		else:
			output = "#2: %.2f\t%.2f\t%d/%d\t%d: %d/%d\t%d: %d/%d\t%s" %(zscore, baseOver_mean, baseOver, self.genomeSize, self.dataA_size, self.getOverA(), self.dataA_nb, self.dataB_size,self.getOverB(),  self.dataB_nb, fileB_name)
		return [output, baseOver_mean, baseOver_l, baseOver, 'BasesOver.', zscore, baseOver_min, baseOver_max]

	########################################################################
	### Calculate the region overlap statistics and return a output_line ###
	########################################################################
	def encodeRegionOver(self, randMatrices, fileB_name, args):
		def calcRegionOver(pAB, pA):
			return pAB / pA

		regionOver = calcRegionOver(self.getOverA(), self.dataA_nb)
		regionOver_l = []
		for m in randMatrices:
			regionOver_l.append( calcRegionOver( m.getOverA(), m.dataA_nb) )
		regionOver_mean = sum(regionOver_l) / len(regionOver_l)
		zscore = self.getZScore(regionOver, regionOver_l)[0]
		regionOver_min = calcRegionOver(0, self.dataA_nb)
		regionOver_max = calcRegionOver(self.dataA_nb, self.dataA_nb)

		if args['-l'] == 'all':
			output = "\tRegions\t%g\t%g\t%g" %(regionOver, zscore, regionOver_mean)
		else:
			output = "#2: %.2f\t%.2f\t%d/%d\t%d: %d/%d\t%d: %d/%d\t%s" %(zscore, regionOver_mean, regionOver, self.genomeSize, self.dataA_size, self.getOverA(), self.dataA_nb, self.dataB_size,self.getOverB(),  self.dataB_nb, fileB_name)
		return [output, regionOver_mean, regionOver_l, regionOver, 'RegionsOver.', zscore, regionOver_min, regionOver_max]


	##################################################################
	### Calculate the pairwise enrichment and return a output_line ###
	##################################################################
	def ernstKellis(self, randMatrices, fileB_name, args):
		def calcErnstKellis(pAB, pA, pB, genomeSize):
			const = 0 # pseudocount of 100 bases for smoothing in the formula
			return (pAB + const) / ((pA * pB)/genomeSize + const)

		pwEnrich = calcErnstKellis(self.overlapSize, self.dataA_size, self.dataB_size, self.genomeSize)
		pwEnrich_l = []
		for m in randMatrices:
			pwEnrich_l.append( calcErnstKellis(m.overlapSize, m.dataA_size, m.dataB_size, m.genomeSize) )
		pwEnrich_mean = sum(pwEnrich_l) / len(pwEnrich_l)
		zscore = self.getZScore(pwEnrich, pwEnrich_l)[0]
		pwEnrich_min = calcErnstKellis(0, self.dataA_size, self.dataB_size, self.genomeSize)
		pwEnrich_max = calcErnstKellis(min(self.dataA_size, self.dataB_size), self.dataA_size, self.dataB_size, self.genomeSize)
		
		if args['-l'] == 'all':
			output = "\tPW Enri\t%g\t%g\t%g" %(pwEnrich, zscore, pwEnrich_mean)
		else:
			output = "#2: %.2f\t%.2f\t%d/%d\t%d: %d/%d\t%d: %d/%d\t%s" %(zscore, pwEnrich_mean, pwEnrich, self.genomeSize, self.dataA_size, self.getOverA(), self.dataA_nb, self.dataB_size,self.getOverB(),  self.dataB_nb, fileB_name)
		return [output, pwEnrich_mean, pwEnrich_l, pwEnrich, 'PairwiseEnri.', zscore, pwEnrich_min, pwEnrich_max]


	#############################################################################
	### Calculate the Pearson correlation coefficent and return a output_line ###
	#############################################################################
	def pearson(self, randMatrices, fileB_name, args):
		# Formula of the coefficient
		def calcPearson(pA, pB, pA_sq, pB_sq, pAB, n):
			numerator = pAB - ((pA * pB) / n)
			denominator = math.sqrt( (pA_sq - pow(pA,2) / n) * (pB_sq - pow(pB,2) / n) )
			if denominator == 0: 
				return 1
			else: 
				return numerator / denominator		

		# Calculte the coefficient
		pearson = calcPearson( self.pA, self.pB, self.pA_sq, self.pB_sq, self.overlapSize, self.genomeSize )
		pearson_l = []
		for m in randMatrices:
				res = calcPearson( m.pA, m.pB, m.pA_sq, m.pB_sq, m.overlapSize, m.genomeSize )
				pearson_l.append( res )
		pearson_mean = sum(pearson_l) / len(pearson_l)
		zscore = self.getZScore(pearson, pearson_l)[0]
		pearson_min = calcPearson( self.pA, self.pB, self.pA_sq, self.pB_sq, 0, self.genomeSize )
		pearson_max = calcPearson( self.pA, self.pB, self.pA_sq, self.pB_sq, min(self.dataA_size, self.dataB_size), self.genomeSize )		

		# Create an output result
		if args['-l'] == 'all':
			output = "\tPearson\t%g\t%g\t%g" %(pearson, zscore, pearson_mean)
		else:
			output = "#2: %.2f\t%.2f\t%d/%d\t%d: %d/%d\t%d: %d/%d\t%s" %(zscore, pearson_mean, pearson, self.genomeSize, self.dataA_size, self.getOverA(), self.dataA_nb, self.dataB_size,self.getOverB(),  self.dataB_nb, fileB_name)
		return [output, pearson_mean, pearson_l, pearson, 'Pearson', zscore, pearson_min, pearson_max]


	######################################################################################
	### Calculate the Normalized Pointwise Mutual Information and return a output_line ###
	######################################################################################
	def npmi(self, randMatrices, fileB_name, args):
		def calcNpmi(pAB, pA, pB, genomeSize):
			pAB = pAB / genomeSize
			pA = pA / genomeSize
			pB = pB / genomeSize
			if pAB == 0 or math.log(pAB) == 0:
				return -1
			else:
				return math.log(pAB/(pA*pB)) / (-math.log(pAB))

		npmi = calcNpmi(self.overlapSize, self.dataA_size, self.dataB_size, self.genomeSize)
		npmi_l = []
		for m in randMatrices:
			npmi_l.append( calcNpmi(m.overlapSize, m.dataA_size, m.dataB_size, m.genomeSize) )
		npmi_mean = sum(npmi_l) / len(npmi_l)
		zscore = self.getZScore(npmi, npmi_l)[0]
		npmi_min = calcNpmi(0, self.dataA_size, self.dataB_size, self.genomeSize)
		npmi_max = calcNpmi(min(self.dataA_size, self.dataB_size), self.dataA_size, self.dataB_size, self.genomeSize)

		if args['-l'] == 'all':
			output = "\tNPMI\t%g\t%g\t%g" %(npmi, zscore, npmi_mean)
		else:
			output = "#2: %.2f\t%.2f\t%d/%d\t%d: %d/%d\t%d: %d/%d\t%s" %(zscore, npmi_mean, npmi, self.genomeSize, self.dataA_size, self.getOverA(), self.dataA_nb, self.dataB_size,self.getOverB(),  self.dataB_nb, fileB_name)
		return [output, npmi_mean, npmi_l, npmi, 'NPMI', zscore, npmi_min, npmi_max]


	######################################################################################
	### Calculate the p-value of binominal test and return a output_line ###
	######################################################################################
	def pBinom(self, randMatrices, fileB_name, args):
		pbinomTest = ro.r['pbinom']
		pbinom = -float( pbinomTest( self.overlapsRegionNumber-1, self.dataA_nb, self.dataB_size/self.genomeSize, lower_tail=args['-f'], log_p=True ).r_repr() ) / math.log(10)
		pbinom_l = []
		for m in randMatrices:
			pbinom_l.append( -float( pbinomTest( m.overlapsRegionNumber, m.dataA_nb, m.dataB_size/m.genomeSize, lower_tail=args['-f'], log_p=True ).r_repr() ) / math.log(10) )
		pbinom_mean = sum(pbinom_l) / len(pbinom_l)
		zscore = 0
		pbinom_min = -float( pbinomTest( 0, self.dataA_nb, self.dataB_size/self.genomeSize, lower_tail=args['-f'], log_p=True ).r_repr() ) / math.log(10)
		pbinom_max = -float( pbinomTest( self.dataA_nb-1, self.dataA_nb, self.dataB_size/self.genomeSize, lower_tail=args['-f'], log_p=True ).r_repr() ) / math.log(10)
		if args['-l'] == 'all':
			output = "\tpBinom\t%g\t%g\t%g" %(pbinom, zscore, pbinom_mean)
		else:
			output = "#2: %.2f\t%.2f\t%d/%d\t%d: %d/%d\t%d: %d/%d\t%s" %(zscore, pbinom_mean, pbinom, self.genomeSize, self.dataA_size, self.getOverA(), self.dataA_nb, self.dataB_size,self.getOverB(),  self.dataB_nb, fileB_name)
		return [output, pbinom_mean, pbinom_l, pbinom, 'pBinom', zscore, pbinom_min, pbinom_max]

