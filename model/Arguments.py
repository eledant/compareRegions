"""compareRegions
-Compare '.bed' or '.bedGraph' files of regions & calculate p-values using randomized data

Usage: compareRegions.py [options] <genome_file> <A_file> <B_files>...

Options:
	-m <arg>	  Number of randomizations for the <A_file> 		[default: 10]
	-n <arg>	  Number of randomizations for each <B_files>  		[default: 10]
	-r 		  Random seed 						[default: False]
	-i <arg>	  Ignore region scores (A|B|AB) 			[default: None]
	-v <arg>	  Verbose output (all|refG|randG|remap|fileA|fileB)	[default: None]
	-l <arg>	  Model (def|all|jac|ebo|ero|pwe|psn|npm)		[default: def]
	-s		  Sort the output by score				[default: False]
	-p		  Plot the results in a graph				[default: False]
	-c		  Create a summary file of the comparaison		[default: False]
	-t <arg>	  Set the number of threads				[default: 1]
"""
from docopt import docopt
import re
import os.path

# 'docopt.py' need to be in the same repository or this option parser need to be install on your computer
class Arguments(dict):

	# Initialize the dictionary with the command line arguments (keys: <genome_file>, <A_file>, <B_files> or "option name")
	def __init__(self):
		self.update( docopt(__doc__, version='1.0') )
		self.checkArguments()

	# Check if the given arguments are syntactically correct
	def checkArguments(self):
		errorList = []
		# Options test
		if not self['-m'].isdigit():
			errorList.append( "'-m' option only takes an integer as argument." )
		elif self['-m'] <= 0:
			errorList.append( "'-m' option only takes an positive integer as argument." )
		if not self['-n'].isdigit():
			errorList.append( "'-n' option only takes an integer as argument." )
		elif self['-n'] <= 0:
			errorList.append( "'-n' option only takes an positive integer as argument." )
		if self['-i'] not in ['A','B','AB','None']:
			errorList.append( "'-i' option only takes 'A', 'B' or 'AB' as argument." )
		if self['-l'] not in ['def','all','jac','ebo','ero','pwe','psn', 'npm']:
			errorList.append( "'-l' option only takes 'def', jac', 'enc', 'pwe', 'psn' or 'npm' as argument." )
		if self['-v'] not in ['None', 'all','refG','randG','remap', 'fileA','fileB']:
			errorList.append( "'-v' option only takes 'all', 'refG', 'randG', 'remap', 'fileA' or 'fileB' as argument." )		

		#Files test
		if not os.path.isfile(self['<genome_file>']):
			errorList.append( "The genome file \"" + self['<genome_file>'] + "\" don't exist or isn't a file." )
		elif re.search('\.bed$', self['<genome_file>']) is None and re.search('\.bedGraph$', self['<genome_file>']) is None:
			errorList.append( "The genome file \"" + self['<genome_file>'] + "\" must be a '.bed' or '.bedGraph' file." )
		if not os.path.isfile(self['<A_file>']):
			errorList.append( "The A file \"" + self['<A_file>'] + "\" don't exist or isn't a file." )
		elif re.search('\.bed$', self['<A_file>']) is None and re.search('\.bedGraph$', self['<A_file>']) is None:
			errorList.append( "The A file \"" + self['<A_file>'] + "\" must be a '.bed' or '.bedGraph' file." )
		for fileB in self['<B_files>']:
			if not os.path.isfile(fileB):
				errorList.append( "The B file \"" + fileB + "\" don't exist or isn't a file." )
			elif re.search('\.bed$', fileB) is None and re.search('\.bedGraph$', fileB) is None:
				errorList.append( "The B file \"" + fileB + "\" must be a '.bed' or '.bedGraph' files." )
		# Print the errors	
		if errorList:
			for e in errorList: print e
			exit()

