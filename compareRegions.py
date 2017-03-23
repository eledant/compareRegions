#!/usr/bin/python
from model import Dataset
from input import Arguments
import pprint, math, datetime


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
		k += 1
		M += (e - tempM)/k
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
	return [ z_score, pval, mean, sd]
	

##################################
### Print output final results ###
##################################
def output_stats(stats, fileB_name):
	print "#3_subject\t", abs(stats[0]), "\t", stats[0], "\t", stats[1], "\t", fileB_name

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
	if args['-v']: refGenome.printReferenceeGenome()

	# Create a Dataset object based on the <A_file>
	fileA = Dataset( args['<A_file>'] )
	fileA.calcDataCoords(refGenome)

	# Print output header
	#output_header(args, fileA)

	# Create n randomizations of the <A_file>
	randFileA = fileA.randomize(args, refGenome)

	# Foreach <B_files>
	for fileB_name in args['<B_files>']:

		# Create a Dataset object based on the <B_file>
		fileB = Dataset( fileB_name )
		fileB.calcDataCoords(refGenome)

		# Compare <A_file> and <B_file>
		res = fileA.compareData(fileB, args)
		overlapBP = res[2]
		exit()

		# Create n randomizations of the <B_file>
		randFileB = fileB.randomize(args, refGenome)

		# Compare randomized <A_file> and randomized <B_file>
		randOverlapBP = []
		for nbRandom in range( int(args['-n']) ):
			res = randFileA.compareData(randFileB, args)
			randOverlapBP.append( res[2] )
			print res

		# Calculate stats
		stats = getStats(overlapBP, randOverlapBP)

		# Print output
		output_stats(stats, fileB_name)



