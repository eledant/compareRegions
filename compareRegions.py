#!/usr/bin/python
from model import Dataset
from input import Arguments
import pprint


# Print output header with columns names
def output_header(args, fileA):
	# Count regions & bp:
	totalRegions, totalBP = 0, 0
	for region in fileA:		
		for data in fileA[region]:
			totalRegions += 1
			totalBP += data['chromEnd'] - data['chromStart'] + 1
	# Output
	print "#1_query\tA_filename=\"",args['<A_file>'],"\"\tA_bp=",totalBP,"\tA_regions=",totalRegions
	print "#2_fields\tabs(z-score)\tz-score\tp-val\tB_name\tB_filename\tbp_overlap(obs/exp)\tB_bp\tB_regions\tA_region_overlaps(obs/exp)\tB_region_overlaps(obs/exp)\tregion_overlap_z-score\tregion_overlap_p-val"

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
		print fileA.compareData(fileB, args)

		# Create n randomizations of the <B_file>
		#randFileB = fileB.randomize(args, refGenome)

		# Compare randomized <A_file> and randomized <B_file>
		#randFileA.compareData(randFileB, args)

		# Calculate stats

		# Print output




