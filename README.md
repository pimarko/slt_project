## Segmentation of Neuron Bundles from Diffusion MRI


### Python, Libraries and install.sh
In order to get all libraries run the following code in terminal/command line

./install.sh


### Module find_clusters.py

At the beginning a library find_clusters is imported. There, the function to get the
clusters is defined.

import find_clusters

### modify constants

This section explains the constants which can be modified.

**K_NN = 26**								
this is number of k nearest neighbors of a spin

**K_NN_SUBSET = 6**									
from these K_NN a subset can be chosen, s.t. the maximal signal strength so only the K_NN_SUBSET neighbors with  highest signal strength are chosen as  neighbors

**VOXELS_GRID = [10,10,10]** 						
define the size of the cube of voxels considerartion of a subset of the brain it also sais from where to where we are going

**Q = 10**
Q is the number of spin states a spin can take as explained in the Blatt 1997 paper

**M = 100**
number of monte carlo iterations exponential decay rate

**GENERATE_PLOT_SEARCH_SPM = True**
set this boolean to True if the search for superparamagnetic phase should be plotted	

**GENERATE_PLOT_CLUSTERS = True**
set this boolean to True if the final clusters should be plotted

**SUBSET_KNN = True**
decide if you want to use K_NN_SUBSET
True: 	use of K_NN_SUBSET as neighbors with the most similiar diffusion profile
False: 	use of K_NN as neighbors

**SMOOTHING = False**
a voxel is assigned to a cluster if its G_ij > 0.5, to get a smoother cluster the voxels at the bounderies are also assigned to the cluster 

**thres = 1**
set a threshold for a number of data point in the cluster which will be plotted > thres take it

**thres_pic = 1**								
threshold for gif

**percent_in_cluster = 5**
The minimal % of elements in clusters that should be fullfiled while finding the optimal T_final

**DISIMILARITY_MEASURE_INNER_PRODUCT = False**
set it to True if the inner product of the signal strenght of two voxels should be considered.

**NORMALIZE_0_1 = False**
this normalizes the data_matrix to 0-1 range, where each signal strength direction is normalized

**Grid Start Position** 
With the next three constants a place of the cube can be defined

init_i = 135, choose at which voxel in i-direction to begin with 
init_j = 76, choose at which voxel in j-direction to begin with
init_k = 74, choose at which voxel in k-direction to begin with







