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
		self.dataA_nb, self.dataB_nb = nA, nB
		# Total size of all the regions for each dataset 
		self.dataA_size, self.dataB_size = sizeA, sizeB


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


	##############################################################################
	### Return the overlaps size and the number of overlap in A and B datasets ###
	##############################################################################
	def countOverlaps(self):
		overA, overB, overBP = 0, 0, 0
		for i in range(len(self)):
			overBP += self[i]
		overA = len(list(set(self.indexA_list))) - len(list(set(self.splitRegionA)))
		overB = len(list(set(self.indexB_list))) - len(list(set(self.splitRegionB)))
		return overA, overB, overBP

	############################################
	### Set the size of the reference genome ###
	############################################
	def addGenomeSize(self, _genomeSize):
		self.genomeSize = _genomeSize

	##########################################################
	### Return the sum of the product bewteen each overlap ###
	##########################################################
	def sumProduct(self, matrice):
		mX, mY = self, matrice
		sumP = 0
		indexY = 0
		for x in range(len(mX)):
			y = indexY
			while y < len(mY):
				if mX.indexA_list[x] == mY.indexA_list[y] and mX.indexB_list[x] == mY.indexB_list[y]:
					sumP += mX[x] * mY[y]
				elif mY.indexA_list[y] > mX.indexA_list[x] or mY.indexB_list[y] > mX.indexB_list[x]:
					indexY = y
					break
				y += 1					
		return sumP

# -----------------------------------------------------------------------------------------

	########################################################################################
	### Menu : select the method to calculate the score according to the selected option ###
	########################################################################################
	def calcStats(self, randMatrices, fileB_name, args):
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
		overlapA, overlapB, overlapBP = self.countOverlaps()
		overlapAB = overlapA * overlapB
		for randMatrix in randMatrices:
			randRes = randMatrix.countOverlaps()
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
		inter = self.countOverlaps()[2]
		union = self.dataA_size + self.dataB_size - inter
		jaccard = float(inter) / float(union)
		
		randJacc = []
		for matrix in randMatrices:
			inter = matrix.countOverlaps()[2]
			union = matrix.dataA_size + matrix.dataB_size - inter
			randJacc.append( float(inter) / float(union) )
		mean = sum(randJacc) / len(randJacc)

		fileB_name = fileB_name.split('/')[-1]
		output = "#3_subject\t%.4e\t%.4e\t%s" %(mean, jaccard, fileB_name)
		return [output, mean]


	#####################################################################################
	### Calculate the region and basepair overlap statistics and return a output_line ###
	#####################################################################################
	def encode(self, randMatrices, fileB_name):
		overA, overB, interBP = self.countOverlaps()
		baseOver = float(interBP) / float(self.dataA_size)
		regionOver = float(overA) / float(self.dataA_nb)

		randBaseOver, randRegionOver = [], []
		for matrix in randMatrices:
			overA, overB, interBP = matrix.countOverlaps()
			randBaseOver.append( float(interBP) / float(matrix.dataA_size) )
			randRegionOver.append( float(overA) / float(matrix.dataA_nb) )
		randBaseOver_mean = sum(randBaseOver) / len(randBaseOver)
		randRegionOver_mean = sum(randRegionOver) / len(randRegionOver)

		fileB_name = fileB_name.split('/')[-1]
		output = "#3_subject\t%.4e\t%.4e\t%s\t%.4e\t%.4e" %(randBaseOver_mean, randRegionOver_mean, fileB_name, baseOver, regionOver)
		return [output, randBaseOver_mean]


	##################################################################
	### Calculate the pairwise enrichment and return a output_line ###
	##################################################################
	def ernstKellis(self, randMatrices, fileB_name):
		overBP = self.countOverlaps()[-1]
		overAB_chance = float(self.dataB_size) / float(self.genomeSize)
		overBP_expct = float(self.dataB_size) / overAB_chance
		pw_enrich = overBP / overBP_expct

		pw_randEnrich = []
		for matrix in randMatrices:
			overBP = matrix.countOverlaps()[-1]
			overAB_chance = float(matrix.dataB_size) / float(matrix.genomeSize)
			overBP_expct = float(matrix.dataB_size) / overAB_chance
			pw_randEnrich.append( overBP / overBP_expct )
		pw_randEnrich_mean = sum(pw_randEnrich) / len(pw_randEnrich)
		
		fileB_name = fileB_name.split('/')[-1]
		output = "#3_subject\t%.4e\t%.4e\t%s" %(pw_randEnrich_mean, pw_enrich, fileB_name)
		return [output, pw_randEnrich_mean]


	#############################################################################
	### Calculate the Pearson correlation coefficent and return a output_line ###
	#############################################################################
	def pearson(self, randMatrices, fileB_name):
		# n = number of comparaison, sumX = sum of X overlap, sumXSq = sum of the Y squared overlap
		n = float(self.dataA_nb * self.dataB_nb)
		sumX = self.countOverlaps()[-1]
		sumXSq = float(sum([pow(self[i],2) for i in range(len(self))]))
		# Do the same thing with a list of matrice + calculate the mean for each stats
		sumY_l, sumYSq_l, pSum_l = [], [], []
		for matrix in randMatrices:
			sumY_l.append( matrix.countOverlaps()[-1] )
			sumYSq_l.append( float(sum( [pow(matrix[i],2) for i in range(len(matrix))])) )
			pSum_l.append( self.sumProduct(matrix) )
		sumY = sum(sumY_l) / len(sumY_l)
		sumYSq = sum(sumYSq_l) / len(sumYSq_l)
		pSum = sum(pSum_l) / len(pSum_l)
		# Calculte the coefficient
		numerator = pSum - ((sumX * sumY) / n)
		denominator = math.sqrt( (sumXSq - pow(sumX,2) / n) * (sumYSq - pow(sumY,2) / n) )
		#print "n =", n
		#print "sumX =", sumX
		#print "sumXSq =", sumXSq			
		#print "sumY =", sumY
		#print "sumYSq =", sumYSq		
		#print "pSum =", pSum
		#print "numerator =", numerator
		#print "denominator =", denominator
		if denominator == 0: 
			pearson = 1
		else: 
			pearson = numerator / denominator			
		# Create an output result
		fileB_name = fileB_name.split('/')[-1]
		output = "#3_subject\t%.2f\t%s" %(pearson, fileB_name)
		return [output, pearson]


	def npmi(self, randMatrices, fileB_name):
		print "TO DO"




