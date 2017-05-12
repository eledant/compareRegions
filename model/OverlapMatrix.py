import math

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


	################################
	### Return the overlaps size ###
	################################
	def getOverlapsSize(self):
		total = 0
		for i in range(len(self)):
			total += self[i]
		return total


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
	
		print "Overlap Size =", self.getOverlapsSize()

		if args['-l'] in ['all', 'def']:
			res_stats.append( self.zScore(randMatrices, fileB_name, args) )
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
		return res_stats

	
	######################################################
	### Calculate the z-score and return a output_line ###
	######################################################
	def zScore(self, randMatrices, fileB_name, args):
		# Return the z-score, p-value, mean and sd
		def getZScore(value, distribution):
			gt, SUM, M, S, k = 0, 0, 0, 0, 1
			for e in distribution:
				if e > value: 
					gt += 1
				SUM += e
				tempM = M
				M += (e - tempM)/k
				k += 1
				S += (e - tempM)*(e - M)
			pval = 2*(gt/(k-1))
			if pval > 1:
				pval = 2 - pval
			mean = SUM/(k-1)
			sd = math.sqrt(S/(k-1))
			if sd != 0:
				z_score = (value - mean) / sd
			else:
				z_score = 'NaN'
			return [z_score, pval, mean, sd]

		# Create output final results
		def output_stats(statsBP, statsAB, fileB_name, args):
			fileB_name = fileB_name.split('/')
			# BP values : 0=z_score ; 1=p_value ; 2=expBP ; 4=overlapBP ; 5=totalBP_B ; 6=totalRegions_B
			p_valueBP = statsBP[1]
			p_strBP = 'p='
			if p_valueBP == 0:
				p_valueBP = 2/math.pow( int(args['-n']), 2 )
				p_strBP = 'p<'
			# AB values : 0=z_score ; 1=p_value ; 2=expAB ; 4=A_exp ; 5=B_exp; 6=overlapA ; 7=overlapB
			p_valueAB = statsAB[1]
			p_strAB = 'p='
			if p_valueAB == 0:
				p_valueAB = 2/math.pow( int(args['-n']), 2 )
				p_strAB = 'p<'	

			output_line = ''

			if args['-l'] == 'all':
				if statsBP[0] == 'NaN':
					output_line = "\tZ-Score\tNaN"
				else:
					output_line = "\tZ-Score\t%g" % abs(statsBP[0])
			else:
				if statsBP[0] != 'NaN':
					output_line += "#3_subject\t%.2f\t%.2f" %(abs(statsBP[0]), statsBP[0])
				else:
					output_line += "#3_subject\t" + statsBP[0] + "\t" + statsBP[0]
				output_line += "\t" + "%s%.3g\t" %(p_strBP, p_valueBP) + fileB_name[-1]
				output_line += "\t%.1f/%.1f\t%d\t%d\t%.1f/%.1f\t%.1f/%.1f\t" %(statsBP[4], statsBP[2], statsBP[5], statsBP[6],statsAB[6], statsAB[4], statsAB[7], statsAB[5])
				if statsAB[0] != 'NaN':
					output_line += "%.2f" % statsAB[0]
				else:
					output_line += statsAB[0]
				output_line += "\t" + "%s%.3g" %(p_strAB, p_valueAB)
			return output_line

		# Get the overlaps size and the number of overlap in A and B datasets
		A_exp, B_exp, randOverlapBP, randOverlapAB = 0, 0, [], []
		overlapA, overlapB, overlapBP = self.getOverA(), self.getOverB(), self.getOverlapsSize()
		overlapAB = overlapA * overlapB
		for randMatrix in randMatrices:
			randRes = randMatrix.getOverA(), randMatrix.getOverB(), randMatrix.getOverlapsSize()
			randOverlapAB.append( randRes[0] * randRes[1] )
			randOverlapBP.append( randRes[2] )
			A_exp += randRes[0]
			B_exp += randRes[1]
		# Calculate the z-score and other stats
		A_exp /= float(args['-m']) * float(args['-n'])
		B_exp /= float(args['-m']) * float(args['-n'])
		statsBP = getZScore(overlapBP, randOverlapBP) + [overlapBP, self.dataB_size, self.dataB_nb]
		statsAB = getZScore(overlapAB, randOverlapAB) + [A_exp, B_exp, overlapA, overlapB]
		# Create output for the <B_file>
		output = output_stats(statsBP, statsAB, fileB_name, args)
		return [output, statsBP[0], randOverlapBP, overlapBP, 'Z-Score']


	#################################################################
	### Calculate the Jaccard similarity and return a output_line ###
	#################################################################
	def jaccard(self, randMatrices, fileB_name, args):
		def calcJaccard(pAB, pA, pB):
			if pA + pB - pAB == 0:
				return 0
			else:
				return pAB / (pA + pB - pAB)
		
		jaccard = calcJaccard(self.getOverlapsSize(), self.dataA_size, self.dataB_size)
		jaccard_l = []
		for m in randMatrices:
			jaccard_l.append( calcJaccard(m.getOverlapsSize(), m.dataA_size, m.dataB_size) )
		jaccard_mean = sum(jaccard_l) / len(jaccard_l)

		if args['-l'] == 'all':
			output = "\tJaccard\t%g\t%g" %(jaccard_mean, jaccard)
		else:
			output = "#3_subject\t%g\t%g\t%s" %(jaccard_mean, jaccard, fileB_name)
		return [output, jaccard_mean, jaccard_l, jaccard, 'Jaccard Similarity']


	##########################################################################
	### Calculate the basepair overlap statistics and return a output_line ###
	##########################################################################
	def encodeBaseOver(self, randMatrices, fileB_name, args):
		def calcBaseOver(pAB, pA):
			return pAB / pA

		baseOver = calcBaseOver(self.getOverlapsSize(), self.dataA_size)
		baseOver_l = []
		for m in randMatrices:
			baseOver_l.append( calcBaseOver( m.getOverlapsSize(), m.dataA_size) )
		baseOver_mean = sum(baseOver_l) / len(baseOver_l)

		if args['-l'] == 'all':
			output = "\tBases\t%g\t%g" %(baseOver_mean, baseOver)
		else:
			output = "#3_subject\t%g\t%g\t%s" %(baseOver_mean, baseOver, fileB_name)
		return [output, baseOver_mean, baseOver_l, baseOver, 'Bases Overlap']

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

		if args['-l'] == 'all':
			output = "\tRegions\t%g\t%g" %(regionOver_mean, regionOver)
		else:
			output = "#3_subject\t%g\t%g" %(regionOver_mean, regionOver, fileB_name)
		return [output, regionOver_mean, regionOver_l, regionOver, 'Regions Overlap']


	##################################################################
	### Calculate the pairwise enrichment and return a output_line ###
	##################################################################
	def ernstKellis(self, randMatrices, fileB_name, args):
		def calcErnstKellis(pAB, pA, pB, genomeSize):
			const = 0 # pseudocount of 100 bases for smoothing in the formula
			return (pAB + const) / ((pA * pB)/genomeSize + const)

		pwEnrich = calcErnstKellis(self.getOverlapsSize(), self.dataA_size, self.dataB_size, self.genomeSize)
		pwEnrich_l = []
		for m in randMatrices:
			pwEnrich_l.append( calcErnstKellis(m.getOverlapsSize(), m.dataA_size, m.dataB_size, m.genomeSize) )
		pwEnrich_mean = sum(pwEnrich_l) / len(pwEnrich_l)
		
		if args['-l'] == 'all':
			output = "\tPW Enri\t%g\t%g" %(pwEnrich_mean, pwEnrich)
		else:
			output = "#3_subject\t%g\t%g\t%s" %(pwEnrich_mean, pwEnrich, fileB_name)
		return [output, pwEnrich_mean, pwEnrich_l, pwEnrich, 'Pairwise Enrichment']


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
		pearson = calcPearson( self.pA, self.pB, self.pA_sq, self.pB_sq, self.getOverlapsSize(), self.genomeSize )
		pearson_l = []
		for m in randMatrices:
				res = calcPearson( m.pA, m.pB, m.pA_sq, m.pB_sq, m.getOverlapsSize(), m.genomeSize )
				pearson_l.append( res )
		pearson_mean = sum(pearson_l) / len(pearson_l)
		
		# Create an output result
		if args['-l'] == 'all':
			output = "\tPearson\t%g\t%g" %(pearson_mean, pearson)
		else:
			output = "#3_subject\t%g\t%g\t%s" %(pearson_mean, pearson, fileB_name)
		return [output, pearson_mean, pearson_l, pearson, 'Pearson']


	######################################################################################
	### Calculate the Normalized Pointwise Mutual Information and return a output_line ###
	######################################################################################
	def npmi(self, randMatrices, fileB_name, args):
		def calcNpmi(pAB, pA, pB, genomeSize):
			pAB = pAB / genomeSize
			pA = pA / genomeSize
			pB = pB / genomeSize
			if pAB == 0 or math.log(pAB) == 0:
				return 0
			else:
				return math.log(pAB/(pA*pB)) / (-math.log(pAB))

		npmi = calcNpmi(self.getOverlapsSize(), self.dataA_size, self.dataB_size, self.genomeSize)
		npmi_l = []
		for m in randMatrices:
			npmi_l.append( calcNpmi(m.getOverlapsSize(), m.dataA_size, m.dataB_size, m.genomeSize) )
		npmi_mean = sum(npmi_l) / len(npmi_l)

		if args['-l'] == 'all':
			output = "\tNPMI\t%g\t%g" %(npmi_mean, npmi)
		else:
			output = "#3_subject\t%g\t%g\t%s" %(npmi_mean, npmi, fileB_name)
		return [output, npmi_mean, npmi_l, npmi, 'NPMI']


