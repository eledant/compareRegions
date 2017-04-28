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


	##############################################
	### Add A and B scores for the new overlap ###
	##############################################
	def addScores(self, a, b):
		self.scoreA.append(a)
		self.scoreB.append(b)


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

# -----------------------------------------------------------------------------------------

	########################################################################################
	### Menu : select the method to calculate the score according to the selected option ###
	########################################################################################
	def calcStats(self, randMatrices, fileB_name, args):
		fileB_name = fileB_name.split('/')[-1]
		if args['-l'] == 'def':
			line, key = self.zScore(randMatrices, fileB_name, args)
		elif args['-l'] == 'jac':
			line, key = self.jaccard(randMatrices, fileB_name)
		elif args['-l'] == 'enc':
			line, key = self.encode(randMatrices, fileB_name)
		elif args['-l'] == 'pwe':
			line, key = self.ernstKellis(randMatrices, fileB_name)
		elif args['-l'] == 'psn':
			line, key = self.pearson(randMatrices, fileB_name)
			lineP, keyP = self.zScore(randMatrices, fileB_name, args)
			print "%g %g %s" %(key, keyP, fileB_name)
		elif args['-l'] == 'npm':
			line, key = self.npmi(randMatrices, fileB_name)
		return [line, key]

	
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
		output_line = output_stats(statsBP, statsAB, fileB_name, args)
		return [output_line, statsBP[0]]


	#################################################################
	### Calculate the Jaccard similarity and return a output_line ###
	#################################################################
	def jaccard(self, randMatrices, fileB_name):
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

		output = "#3_subject\t%g\t%g\t%s" %(jaccard_mean, jaccard, fileB_name)
		return [output, jaccard_mean]


	#####################################################################################
	### Calculate the region and basepair overlap statistics and return a output_line ###
	#####################################################################################
	def encode(self, randMatrices, fileB_name):
		def calcEncode(pAB, pA):
			return pAB / pA

		baseOver = calcEncode(self.getOverlapsSize(), self.dataA_size)
		regionOver = calcEncode(self.getOverA(), self.dataA_nb)

		baseOver_l, regionOver_l = [], []
		for m in randMatrices:
			baseOver_l.append( calcEncode( m.getOverlapsSize(), m.dataA_size) )
			regionOver_l.append(  calcEncode( self.getOverA(), m.dataA_nb) )
		baseOver_mean = sum(baseOver_l) / len(baseOver_l)
		regionOver_mean = sum(regionOver_l) / len(regionOver_l)

		output = "#3_subject\t%g\t%g\t%s\t%g\t%g" %(baseOver_mean, regionOver_mean, fileB_name, baseOver, regionOver)
		return [output, baseOver_mean]


	##################################################################
	### Calculate the pairwise enrichment and return a output_line ###
	##################################################################
	def ernstKellis(self, randMatrices, fileB_name):
		def calcErnstKellis(pAB, pA, pB, genomeSize):
			#const = 100 # pseudocount of 100 bases for smoothing
			return (pAB + const) / ((pA * pB)/genomeSize + const)

		pwEnrich = calcErnstKellis(self.getOverlapsSize(), self.dataA_size, self.dataB_size, self.genomeSize)
		pwEnrich_l = []
		for m in randMatrices:
			pwEnrich_l.append( calcErnstKellis(m.getOverlapsSize(), m.dataA_size, m.dataB_size, m.genomeSize) )
		pwEnrich_mean = sum(pwEnrich_l) / len(pwEnrich_l)
		
		output = "#3_subject\t%g\t%g\t%s" %(pwEnrich_mean, pwEnrich, fileB_name)
		return [output, pwEnrich_mean]


	#############################################################################
	### Calculate the Pearson correlation coefficent and return a output_line ###
	#############################################################################
	def pearson(self, randMatrices, fileB_name):
		# Return the sum of a list of number
		def getSum(var):
			size = 0
			for e in var:
				size += e
			return size
		# Return the sum of each squared number
		def getSqrtSum(var):
			size = 0
			for e in var:
				size += pow(e,2)
			return size
		# Return the sum of the multiplication of two list 
		def getMultSum(varA, varB):
			size = 0
			for i in range(len(varA)):
				size += varA[i] * varB[i]
			return size
		# Formula of the coefficient
		def calcPearson(pA, pB, pA_sq, pB_sq, pAB, n):
			numerator = pAB - ((pA * pB) / n)
			denominator = math.sqrt( (pA_sq - pow(pA,2) / n) * (pB_sq - pow(pB,2) / n) )
			if denominator == 0: 
				return 1
			else: 
				return numerator / denominator		

		# Calculte the coefficient
		pearson = calcPearson( getSum(self.scoreA), getSum(self.scoreB), getSqrtSum(self.scoreA), getSqrtSum(self.scoreB), getMultSum(self.scoreA, self.scoreB), len(self.scoreA) )
		pearson_l = []
		for m in randMatrices:
				res = calcPearson( getSum(m.scoreA), getSum(m.scoreB), getSqrtSum(m.scoreA), getSqrtSum(m.scoreB), getMultSum(m.scoreA, m.scoreB), len(m.scoreA) )
				pearson_l.append( res )
		pearson_mean = sum(pearson_l) / len(pearson_l)
		
		# Create an output result
		output = "#3_subject\t%g\t%g\t%s" %(pearson_mean, pearson, fileB_name)
		return [output, pearson]


	######################################################################################
	### Calculate the Normalized Pointwise Mutual Information and return a output_line ###
	######################################################################################
	def npmi(self, randMatrices, fileB_name):
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

		output = "#3_subject\t%g\t%g\t%s" %(npmi_mean, npmi, fileB_name)
		return [output, npmi_mean]
	


