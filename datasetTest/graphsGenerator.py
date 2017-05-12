#!/usr/bin/python
import matplotlib.pyplot as plt
import numpy as np 


zscore = np.array([0, 0, 8.6392, 9.06308, 28.2167, 37.1361, 52.88])
pearson = np.array([-0.000075031, -0.000075031, 0.000040427, 0.00992572, 0.0000260766, 0.0000359273, 0.0000250765])
npmi = np.array([0, 0, 0.00360779, 0.0200667, 0.0150014, 0.0213235, 0.0232104])
jaccard = np.array([0, 0, 0.0000580603, 0.000097908, 0.000505823, 0.0000554952, 0.00005])
bases = np.array([0, 0, 0.00011545, 0.0001952, 0.0001011, 0.00011095, 0.0001001])
regions = np.array([0, 0, 0.0002, 0.0004, 0.00022, 0.00023, 0.000195])
pairwise = np.array([0, 0, 1.53881, 2.60179, 1.34754, 1.47883, 1.33422])

suptitle = '1% overlap | different regions number | 200 000 bp'
scale = np.array([1, 10, 50, 100, 500, 1000, 2000])
allStats = [zscore, pearson, npmi, jaccard, bases, regions, pairwise]
titles = ['Z-Score', 'Pearson', 'NPMI', 'Jaccard', 'Bases Overlap', 'Regions Overlap', 'Pairwise Enrichment']

plt.figure(1)
plt.suptitle(suptitle, fontsize=14, fontweight='bold')
for i in range(len(allStats)):
	plt.subplot(3,3,i+1)
	plt.title(titles[i])
	plt.plot(scale, allStats[i], marker='o')
	plt.xlim([-100, scale[-1]])
manager = plt.get_current_fig_manager()
manager.resize(*manager.window.maxsize())
plt.subplots_adjust(left=0.1, right=0.9, top=0.9, bottom=0.1)
plt.show()
