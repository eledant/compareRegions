# compareRegions
Compare '.bed' or '.bedGraph' files of regions & calculate p-values using randomized data

Usage: compareRegions.py [options] <genome_file> <A_file> <B_files>...

Options:
	-n <arg>	  Number of randomizations of each dataset to analyze 	[default: 10]
	-r 		  Random seed 						[default: False]
	-i <arg>	  Ignore region scores (A|B|AB) 			[default: None]
	-v <arg>	  Verbose output (all|refG|randG|remap|fileA|fileB)	[default: None]
	-m <arg>	  Model (def)						[default: def]
