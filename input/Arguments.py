"""compareRegions
-Compare '.bed' or '.bedGraph' files of regions & calculate p-values using randomized data

Usage: compareRegions.py [options] <genome_file> <A_file> <B_files>...

Options:
	-n <arg>	  Number of randomizations of each dataset to analyze 	[default: 10]
	-r 		  Random seed 						[default: False]
	-i <arg>	  Ignore region scores (A|B|AB) 			[default: None]
	-v 		  Verbose output 					[default: False]
	-m <arg>	  Model (nil|rel|rel_reg|pos|pos_sample)		[default: rel]
	-s <arg>	  Sample region size (required for -m pos_sample) 	[default: 1000]

"""
from docopt import docopt
import re
import os.path

# 'docopt.py' need to be in the same repository or this option parser need to be install on your computer
class Arguments(dict):

	# Initialize the dictionary with the command line arguments (keys: <genome_file>, <A_file>, <B_files> or "option name")
	def __init__(self):
		self.update( docopt(__doc__, version='1.8') )
		self.checkArguments()

	# Check if the given arguments are syntactically correct
	def checkArguments(self):
		errorList = []
		# Options test
		if not self['-n'].isdigit():
			errorList.append( "'-n' option only takes an integer as argument." )
		elif self['-n'] <= 0:
			errorList.append( "'-n' option only takes an positive integer as argument." )
		if self['-i'] not in ['A','B','AB','None']:
			errorList.append( "'-i' option only takes 'A', 'B' or 'AB' as argument." )
		if self['-m'] not in ['nil','rel','rel_reg','pos','pos_sample']:
			errorList.append( "'-m' option only takes 'nil', 'rel', 'rel_reg', 'pos' or 'pos_sample' as argument." )
		elif self['-m'] == 'pos_sample':
			if not self['-s'].isdigit():
				errorList.append( "'-s' option only takes an integer as argument." )
		#Files test
		if not os.path.isfile(self['<genome_file>']):
			errorList.append( "<genome_file> don't exist or isn't a file." )
		elif re.search('\.bed$', self['<genome_file>']) is None and re.search('\.bedGraph$', self['<genome_file>']) is None:
			errorList.append( "<genome_file> must be a '.bed' or '.bedGraph' file." )
		if not os.path.isfile(self['<A_file>']):
			errorList.append( "<A_file> don't exist or isn't a file." )
		elif re.search('\.bed$', self['<A_file>']) is None and re.search('\.bedGraph$', self['<A_file>']) is None:
			errorList.append( "<A_file> must be a '.bed' or '.bedGraph' file." )
		for fileB in self['<B_files>']:
			if not os.path.isfile(fileB):
				errorList.append( "<B_files> don't exist or isn't a file." )
				break
			elif re.search('\.bed$', fileB) is None and re.search('\.bedGraph$', fileB) is None:
				errorList.append( "<B_files> must be a '.bed' or '.bedGraph' files." )
				break
		# Print the errors	
		if errorList:
			for e in errorList: print e
			exit()

