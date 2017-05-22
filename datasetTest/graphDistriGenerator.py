#!/usr/bin/python
import matplotlib.pyplot as plt
import numpy as np
import os
import seaborn as sns
import pandas as pd


folderPath = './'
folderName = 'distri'

distriAll, scale, = {}, []
titles = ['Score Product', 'Pearson', 'NPMI', 'Jaccard', 'BasesOver.', 'RegionsOver.', 'PairwiseEnri.']
fileNames = os.listdir(folderPath+folderName)
fileIDs = [int(x.split('_')[1]) for x in fileNames]
fileIDs.sort()
fileNames = ['distri_'+str(x) for x in fileIDs]

for fileName in fileNames:
	scale.append( fileName.split('_')[1] )
	FILE = open(folderPath+folderName+'/'+fileName, 'r')
	for line in FILE:
		line = line.rstrip('\n')
		statName, distriLine = line.split('\t')
		distri = distriLine.split(',')
		distri = [float(i) for i in distri]
		if statName not in distriAll:
			distriAll[statName] = [ distri ]
		else:
			distriAll[statName].append( distri )
	FILE.close()



plt.figure(1)
suptitle = '1% overlap | different regions number | 20 000 000 bp'
plt.suptitle(suptitle, fontsize=18, fontweight='bold')
i = 0
cm = plt.get_cmap('gist_rainbow')
for statName in titles:
	plt.subplot(3,3,i+1).set_color_cycle([cm(1.*j/len(distriAll[statName])) for j in range(len(distriAll[statName]))])
	plt.title(statName, fontsize=14, fontweight='bold')
	for x in range(len(distriAll[statName])):
		sns.distplot(distriAll[statName][x], hist=False, rug=False, label=scale[x])
	plt.legend(ncol=2)
	i += 1
	
manager = plt.get_current_fig_manager()
manager.resize(*manager.window.maxsize())
plt.subplots_adjust(left=0.1, right=0.9, top=0.9, bottom=0.1)
plt.show()
