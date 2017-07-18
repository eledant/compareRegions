#!/usr/bin/python
from model import Dataset, Arguments
import matplotlib.pyplot as plt
import seaborn as sns
import pprint, datetime, os.path, math
#print datetime.datetime.now().time()

##############################################
### Print output header with columns names ###
##############################################
def output_header(args, fileA_name):
	if args['-l'] != 'all':
		print "#1: z-score\tmean\toverlap/genome\tsizeA: overlapA/nbA\tsizeB: overlapB/nbB\tdatasetB"
	else:
		print "datasetB\tStatName\tScore\tZ-Score\tMean"


#############################################
### Print output_lines sorted output_keys ###
#############################################
def output_results(total_stats, fileB_name):
	keys, lines = [], []
	# Get the keys (scores) and their output lines
	for stats in total_stats:
		lines.append( stats[0] )
		keys.append( stats[1] )
	# Sort the output lines by highest scores
	if args['-s']:
		for i in range(len(keys)):
			keys[i] = abs(keys[i])
		lines = [x for (y,x) in sorted(zip(keys,lines), reverse=True)]
	# Print lines
	if args['-l'] == 'all':
		print ">>> %s" %fileB_name.split('/')[-1]
	for line in lines:
		print line


######################################
### Plot a graph of a distribution ###
######################################
def plotDistributions(total_stats, args, nbRegionA, nbRegionB, bpTotalA, bpTotalB):
	# Need to set manually the title
	prcOver = '1'
	suptitle = '%s%% overlap | %s/%s regions | %s/%s bp | Random region/overlap size | Test 8' %(prcOver, nbRegionA, nbRegionB, bpTotalA, bpTotalB)
	# Create 1 figure for each stat
	sns.set(color_codes=True)
	plt.figure(1)
	for i in range(len(total_stats)):
		if args['-l'] == 'all':
			plt.subplot(3,3,i+1)
		plt.suptitle(suptitle, fontsize=18, fontweight='bold')
		comment = '%s expected = %g | Observed value = %g' % (total_stats[i][4], total_stats[i][1], total_stats[i][3])
		# Check if matrix is not singular
		if len(set(total_stats[i][2])) > 1:
			sns.distplot(total_stats[i][2], kde=True, rug=False)
		plt.title(comment, fontweight='bold')
		plt.annotate('Z-Score = %g\nMin = %g\nMax = %g' %(total_stats[i][5], total_stats[i][6], total_stats[i][7]), xy=(0.80, 0.80), xycoords='axes fraction', bbox=dict(boxstyle="round", fc="w"))
	manager = plt.get_current_fig_manager()
	manager.resize(*manager.window.maxsize())
	plt.subplots_adjust(left=0.1, right=0.9, top=0.9, bottom=0.1)
	plt.show()


######################################################
### Create a summary file of all the stats results ###
######################################################
def summaryStatsOutput(total_stats, nbRegionA, nbRegionB, bpTotalA, bpTotalB, genomeRegion, genomeSize, overlapSize):
	fileName = 'statsSummary.txt'
	lines = []
	# If the summary stats files already exists, read it
	if os.path.exists(fileName):
		FILE = open(fileName, 'r')
		for line in FILE:
			lines.append( line.rstrip('\n') )
		FILE.close()
	# Rewrite output file with former and new results
	FILE = open(fileName, 'w')
	for i in range(len(total_stats)):
		if total_stats[i][4] == 'Score Product':
			score = total_stats[i][5]
			if lines:
				FILE.write( '%s, %f\n' %(lines[i], score) )
			else:
				FILE.write( 'Z-Score\t%f\n' %score )
		else:
			if lines:
				FILE.write( '%s, %g\n' %(lines[i], total_stats[i][3]) )
			else:
				FILE.write( '%s\t%g\n' %(total_stats[i][4], total_stats[i][3]) )
	FILE.close()


###########################################################
### Create a summary file of all the stats distribution ###
###########################################################
def saveDistribution(total_stats, nbRegionA, nbRegionB, bpTotalA, bpTotalB):
	fileName = 'distri_%s_%s' %(bpTotalA, bpTotalB)
	lines = []
	FILE = open(fileName, 'w')
	for i in range(len(total_stats)):
		FILE.write( '%s\t' %total_stats[i][4] )
		for j in range(len(total_stats[i][2])-1):
			FILE.write( '%f,' %total_stats[i][2][j] )
		FILE.write( '%f\n' %total_stats[i][2][-1] )
	FILE.close()

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
	output_header(args, args['<A_file>'])

	# Foreach <B_files>
	if args['-c']: FH = open('testFiles/summary', 'w')
	mnRandom = 0
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
				res = [randMatrices[-1].getOverA(),randMatrices[-1].getOverB(), randMatrices[-1].overlapSize]
				if args['-c']: FH.write( '%s\t\trandFileA_%d randFileB_%d\t\tbp_overlap: %d\t\tA_region_overlaps: %d\t\tB_region_overlaps: %d\n' \
						%(fileB_name, mRandom, mnRandom, res[2], res[0], res[1]) )
				mnRandom += 1
		 
		# Calculate stats and save results (0=line, 1=score_mean, 2=distribution, 3=expected_score, 4=title)
		total_stats = matrix.calcStats(randMatrices, fileB_name, args)
		# Print output_lines sorted output_keys (z_score or something else)
		output_results(total_stats, fileB_name)
		if args['-l'] == 'all': 
			summaryStatsOutput(total_stats, fileA.getSize()[0], fileB.getSize()[0], fileA.getSize()[1], fileB.getSize()[1], refGenome.getSize()[0], refGenome.getSize()[1], matrix.overlapSize)
			saveDistribution(total_stats, fileA.getSize()[0], fileB.getSize()[0], fileA.getSize()[1], fileB.getSize()[1])
		# Create distribution graph of stats for each file 
		if args['-p']: plotDistributions(total_stats, args, fileA.getSize()[0], fileB.getSize()[0], fileA.getSize()[1], fileB.getSize()[1])

	if args['-c']:  FH.close()


