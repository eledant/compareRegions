#!/usr/bin/python
from model import Dataset, Arguments
import pprint, datetime
#print datetime.datetime.now().time()

##############################################
### Print output header with columns names ###
##############################################
def output_header(args, fileA):
	totalRegions, totalBP = fileA.getSize()
	print "#1_query\tA_filename=\"",args['<A_file>'],"\"\tA_bp=",totalBP,"\tA_regions=",totalRegions
	if args['-l'] == 'def':
		print "#2_fields\tabs(z-score)\tz-score\tp-val\tB_filename\tbp_overlap(obs/exp)\tB_bp\tB_regions\tA_region_overlaps(obs/exp)\tB_region_overlaps(obs/exp)\tregion_overlap_z-score\tregion_overlap_p-val"
	elif args['-l'] == 'jac':
		print "#2_fields\tJaccard similarity\tSimilarity without rand\tB_filename"
	elif args['-l'] == 'enc':
		print "#2_fields\tBase Overlap\tRegion Overlap\tB_filename\tBase Overlap without rand\tRegion Overlap without rand"
	elif args['-l'] == 'pwe':
		print "#2_fields\tPairwise Enrichment\tPairwise Enrichment without rand\tB_filename"
	elif args['-l'] == 'psn':
		print "#2_fields\tPearson Corr Coeff\tPearson Corr Coeff without rand\tB_filename"
	elif args['-l'] == 'npm':
		print "#2_fields\tNPMI\tNPMI without ran\tB_filename"


#############################################
### Print output_lines sorted output_keys ###
#############################################
def output_results(lines, keys):
	for i in range(len(keys)):
		if keys[i] != 'NaN':
			keys[i] = abs(keys[i])
		else:
			keys[i] = 0
	lines = [x for (y,x) in sorted(zip(keys,lines), reverse=True)]
	for line in lines:
		print line

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
	if args['-v'] in ['all', 'refG']: refGenome.printDataset('refGenome')

	# Create a Dataset object based on the <A_file>
	fileA = Dataset( args['<A_file>'] )
	fileA.calcDataCoords(refGenome)
	if args['-v'] in ['all', 'fileA']: fileA.printDataset('fileA')

	# Print output header for the <A_file>
	output_header(args, fileA)

	# Foreach <B_files>
	FH = open('testFiles/summary', 'w')
	output_lines, output_keys, mnRandom = [], [], 0
	for fileB_name in args['<B_files>']:

		# Create a Dataset object based on the <B_file>
		fileB = Dataset( fileB_name )
		fileB.calcDataCoords(refGenome)
		if args['-v'] in ['all', 'fileB']: fileB.printDataset('fileB')

		# Compare <A_file> and <B_file>
		matrix = fileA.compareData(fileB, args)
		matrix.addGenomeSize( refGenome.getSize()[1] )

		randMatrices = []
		# Foreach M number of randomizations
		for mRandom in range( int(args['-m']) ):

			# Create a randomization of the <A_file>
			randFileA = fileA.randomize(refGenome, args, 'randGenomeA_%d' %mRandom)
			if args['-v'] in ['all', 'remap']: randFileA.printDataset('randFileA_%d' %mRandom)

			# Foreach N number of randomizations
			for nRandom in range( int(args['-n']) ):

				# Create a randomization of the <B_file>
				randFileB = fileB.randomize(refGenome, args, 'randGenomeB_%d' %mnRandom)
				if args['-v'] in ['all', 'remap']: randFileB.printDataset('randFileB_%d' %mnRandom)

				# Compare randomized <A_file> and randomized <B_file>
				randMatrices.append( randFileA.compareData(randFileB, args) )
				randMatrices[-1].addGenomeSize( refGenome.getSize()[1] )

				# Print a new line in the summary file
				res = [randMatrices[-1].getOverA(),randMatrices[-1].getOverB(), randMatrices[-1].getOverlapsSize()]
				FH.write( '%s\t\trandFileA_%d randFileB_%d\t\tbp_overlap: %d\t\tA_region_overlaps: %d\t\tB_region_overlaps: %d\n' %(fileB_name, mRandom, mnRandom, res[2], res[0], res[1]) )
				mnRandom += 1
		 
		# Calculate stats
		line, score = matrix.calcStats(randMatrices, fileB_name, args)
		# Create output for the <B_file>
		output_lines.append( line )
		output_keys.append( score )

	FH.close()
	
	# Print output_lines sorted output_keys (z_score or something else)
	#output_results(output_lines, output_keys)


