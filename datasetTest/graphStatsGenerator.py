#!/usr/bin/python
from mpl_toolkits.axes_grid1 import host_subplot
import mpl_toolkits.axisartist as AA
import matplotlib.pyplot as plt
import numpy as np

filePath = './'
fileName = 'statsSummary.txt'

allStats, scale = {}, []
FILE = open(filePath+fileName, 'r')
for line in FILE:
	line = line.rstrip('\n')
	statName, distriLine = line.split('\t')
	distri = distriLine.split(',')
	distri = [float(i) for i in distri]
	allStats[statName] = distri[:]
FILE.close()

scale = np.array(['10', '25', '50', '100', '150', '200', '250', '500', '750', '1000', '1250', '1500', '1750', '2000', '2500', '3000', '4000', '5000', '6000'])
titles = ['Z-Score', 'Pearson', 'NPMI', 'Jaccard', 'BasesOver.', 'RegionsOver.', 'PairwiseEnri.']

plt.figure(1)
suptitle = '1% overlap | 250 regions in A / 10->6000 regions in B | 20 000 000bp | Distri Method 0 | Test4'
plt.suptitle(suptitle, fontsize=18, fontweight='bold')
i = 0
for statsName in titles:
	plt.subplot(3,3,i+1)
	plt.title(statsName)
	plt.xticks(np.arange(0, 6000, 250))
	plt.plot(scale, allStats[statsName], marker='o')
	i += 1
manager = plt.get_current_fig_manager()
manager.resize(*manager.window.maxsize())
plt.subplots_adjust(left=0.1, right=0.9, top=0.9, bottom=0.1)
plt.show()

# -----------------------------------------------------------------------------
# All stats in one graph
plt.suptitle(suptitle, fontsize=25, fontweight='bold')
host = host_subplot(111, axes_class=AA.Axes)
plt.subplots_adjust(right=0.70)
scaleMin = 0
scaleMax = int(scale[-1])*1.2

host.set_xlabel("Number of Region in DatasetB")
parts= []
offset = 80
i = 1
firstPart = True
for statsName in titles:
	if statsName == 'Z-Score':
		host.set_ylabel(statsName)
		p, = host.plot(scale, allStats[statsName], label=statsName, linewidth=6)
		host.axis["left"].label.set_color(p.get_color())
		host.set_xlim(scaleMin, scaleMax)
		host.axis["left"].label.set_size(25)
		host.axis["bottom"].label.set_size(25)
	else:
		part = host.twinx()
		if firstPart:
			firstPart = False
		else:
			new_fixed_axis = part.get_grid_helper().new_fixed_axis
			part.axis["right"] = new_fixed_axis(loc="right", axes=part,offset=(offset*i, 0))
			part.axis["right"].toggle(all=True)
			i += 1
		part.set_ylabel(statsName)
		p, = part.plot(scale, allStats[statsName], label=statsName, linewidth=6)
		part.axis["right"].label.set_color(p.get_color())
		part.set_xlim(scaleMin, scaleMax)
		part.axis["right"].label.set_size(25)
host.legend()
plt.draw()
plt.show()







