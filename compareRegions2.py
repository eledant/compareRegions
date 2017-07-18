#!/usr/bin/python

"""compareRegions
Compare '.bed' or '.bedGraph' files of regions & calculate p-values using randomized data

Usage: compareRegions.py [options] <genomeFile> <fileA> <filesB>...

Options:
-m <arg>	  Number of randomizations for the <fileA> 		[default: 10]
-n <arg>	  Number of randomizations for each <filesB>  		[default: 10]
"""

from docopt import docopt
import pybedtools as pbt
import numpy as np
import pprint, datetime

###########################
###	 START HERE	###
###########################
if __name__ == '__main__':

	print datetime.datetime.now().time()

	# Get the command line arguments in a dictionary
	args = docopt(__doc__, version='1.0')

	# Import the query dataset A
	fileA = pbt.BedTool( args['<fileA>'] )

	# Foreach <B_files>
	for fileB_name in args['<filesB>']:

		# Import one dataset B
		fileB = pbt.BedTool( fileB_name )

		# Compare fileA with fileB to get the observed value
		res = fileA.intersect(fileB)
		obs = sum(float(e[2])-float(e[1])+1 for e in res)
		exp = []

		# Foreach M number of randomizations
		for mRandom in range( int(args['-m']) ):

			# Create a randomization of fileA
			randFileA = fileA.shuffle(g=args['<genomeFile>'])

			# Foreach N number of randomizations
			for nRandom in range( int(args['-n']) ):

				# Create a randomization of fileB
				randFileB = fileB.shuffle(g=args['<genomeFile>'])
		
				# Compare randFileA with randFileB to get the expected values
				res = randFileA.intersect(randFileB)
				exp.append( sum(float(e[2])-float(e[1])+1 for e in res) )

		# Calculate the Z-Score
		mean = np.mean(np.array( exp ))
		sd = np.std(np.array( exp ))
		if sd != 0:
			zscore = (obs - mean) / sd
			print obs, mean, sd
		else:
			zscore = 0
		print 'Z-Score =', zscore

	print datetime.datetime.now().time()

