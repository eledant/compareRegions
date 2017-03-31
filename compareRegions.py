#!/usr/bin/python
from model import Dataset
from input import Arguments
import pprint, math, datetime
#print datetime.datetime.now().time()

##############################################
### Print output header with columns names ###
##############################################
def output_header(args, fileA):
	# Count regions & bp:
	totalRegions, totalBP = 0, 0
	for region in fileA:		
		for data in fileA[region]:
			totalRegions += 1
			totalBP += data['chromEnd'] - data['chromStart'] + 1
	# Output
	print "#1_query\tA_filename=\"",args['<A_file>'],"\"\tA_bp=",totalBP,"\tA_regions=",totalRegions
	print "#2_fields\tabs(z-score)\tz-score\tp-val\tB_filename\tbp_overlap(obs/exp)\tB_bp\tB_regions\tA_region_overlaps(obs/exp)\tB_region_overlaps(obs/exp)\tregion_overlap_z-score\tregion_overlap_p-val"


#####################################################
### Calculate the z-score, p-value, mean and sd?  ###
#####################################################
def getStats(value, distribution):
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
	

##################################
### Print output final results ###
##################################
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


# -----------------------------------------------------------------------------------------

###########################
###	 START HERE	###
###########################
if __name__ == '__main__':
	# Get the command line arguments in a dictionary
	args = Arguments()

	# Create a Dataset object based on the <genome_file>
	refGenome = Dataset( args['<genome_file>'] )
	refGenome.calcGenomeCoords()
	if args['-v'] in ['all', 'refG']: refGenome.printReferenceGenome()

	# Create a Dataset object based on the <A_file>
	fileA = Dataset( args['<A_file>'] )
	fileA.calcDataCoords(refGenome)
	if args['-v'] in ['all', 'fileA']: fileA.printDataset(refGenome, args['<A_file>'])

	# Print output header for the <A_file>
	output_header(args, fileA)

	# Create n randomizations of the <A_file>
	randFileA = fileA.randomize(args, refGenome)

	# Foreach <B_files>
	output_lines, output_keys = [], []
	for fileB_name in args['<B_files>']:

		# Create a Dataset object based on the <B_file>
		fileB = Dataset( fileB_name )
		fileB.calcDataCoords(refGenome)
		if args['-v'] in ['all', 'fileB']: fileB.printDataset(refGenome, fileB_name)

		# Compare <A_file> and <B_file>
		overlapA, overlapB, overlapBP =  fileA.compareData(fileB, args)
		overlapAB = overlapA * overlapB

		# Foreach number of randomizations
		A_exp, B_exp, randOverlapBP, randOverlapAB = 0, 0, [], []
		for nbRandom in range( int(args['-n']) ):

			# Create n randomizations of the <B_file>
			randFileB = fileB.randomize(args, refGenome)

			# Compare randomized <A_file> and randomized <B_file>
			res = randFileA.compareData(randFileB, args)
			randOverlapAB.append( res[0] * res[1] )
			randOverlapBP.append( res[2] )
			A_exp += res[0]
			B_exp += res[1]
		
		# Calculate stats
		totalRegions_B, totalBP_B = 0, 0
		for region in fileB:		
			for data in fileB[region]:
				totalRegions_B += 1
				totalBP_B += data['chromEnd'] - data['chromStart'] + 1
		A_exp /= int(args['-n'])
		B_exp /= int(args['-n'])
		statsBP = getStats(overlapBP, randOverlapBP) + [overlapBP, totalBP_B, totalRegions_B]
		statsAB = getStats(overlapAB, randOverlapAB) + [A_exp, B_exp, overlapA, overlapB]

		# Create output for the <B_file>
		output_lines.append( output_stats(statsBP, statsAB, fileB_name, args) )
		output_keys.append(statsBP[0])

	# Print output_lines sorted output_keys (z_score)
	for i in range(len(output_keys)):
		if output_keys[i] != 'NaN':
			output_keys[i] = abs(output_keys[i])
		else:
			output_keys[i] = 0
	output_lines = [x for (y,x) in sorted(zip(output_keys,output_lines), reverse=True)]
	for line in output_lines:
		print line


