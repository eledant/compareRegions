#!/usr/bin/python
from mpl_toolkits.axes_grid1 import host_subplot
import mpl_toolkits.axisartist as AA
import matplotlib.pyplot as plt
import numpy as np
import math

# Main Parameters
filePath = './'
fileName = 'statsSummary.txt'
titles = ['Z-Score', 'Pearson', 'pBinom', 'NPMI', 'PairwiseEnri.', 'Jaccard', 'BasesOver.', 'RegionsOver.']
colors = ['b', 'r', 'g', 'c', 'm', 'y', 'k', '0.75']

# Label parameters
suptitle = '100% overlap | 1 region of the cluster size | Fix mean score (500 000bp) | Standard comparaison'
xLabel = 'Overlap Percentage'

# Scale paarameters
scale = np.array(['0', '0.5', '1.0', '2.5', '5.0', '7.5', '10.0', '20.0', '30.0', '40.0', '50.0', '75.0', '100.0'])
xScaleStart = 0
xEndScale = 100.1
ticks = 10

# Open file and save measures
allStats = {}
FILE = open(filePath+fileName, 'r')
for line in FILE:
	line = line.rstrip('\n')
	statName, distriLine = line.split('\t')
	distri = distriLine.split(',')
	distri = [float(i) for i in distri]
	allStats[statName] = distri[:]
FILE.close()

# Generate graphs
plt.figure(1)
plt.suptitle(suptitle, fontsize=18, fontweight='bold')
i = 0
for statsName in titles:
	plt.subplot(3,3,i+1)
	plt.xticks(np.arange(xScaleStart, xEndScale, ticks))
	plt.plot(scale, allStats[statsName], marker='o', color=colors[i])
	plt.xlabel(xLabel, color=colors[i], fontsize=18)
	plt.ylabel(statsName, fontsize=24, color=colors[i])
	i += 1
manager = plt.get_current_fig_manager()
manager.resize(*manager.window.maxsize())
plt.subplots_adjust(left=0.1, right=0.9, top=0.9, bottom=0.1)
plt.show()







