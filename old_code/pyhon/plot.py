from time import time

# these two lines have to be inserted here
#import matplotlib
#matplotlib.use('Agg')

#import libraries
import matplotlib,numpy
import nibabel as nib
import numpy as np
import matplotlib.pyplot as plt
import sys
from sklearn.neighbors import BallTree
from math import exp
from random import randint
from collections import Counter
from mpl_toolkits.mplot3d import Axes3D
import itertools
import find_clusters
from scipy.spatial import distance
import matplotlib

thres = 1

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
colors = []
for name, hex in matplotlib.colors.cnames.iteritems():
	colors.append(str(name))

i = 0
for key in new_clusters.keys():
	values = new_clusters[key]
	if(len(values) > thres):
		iss = data_matrix[values,0]
		jss = data_matrix[values,1]
		kss = data_matrix[values,2]

		print(colors[i])
		ax.set_title('Temp = ' + str(T) + ', #clusters= ' + str(len(clusters.keys())), y=1.08)
		ax.scatter(iss, jss, kss,c = colors[i], marker = ',',s = 1000)
		i = (i+3)%len(colors)


ax.set_xlabel('i')
ax.set_ylabel('j')
ax.set_zlabel('k')

#plt.savefig('plot_final_clusters.png')
plt.show()