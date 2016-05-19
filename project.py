#import libraries
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
import seaborn as sns

###cij wrong!!

#modify constants
K_NN = 5
VOXELS_GRID = [7,7,7]
Q = 10
M = 100
eta = 0.97
GENERATE_PLOT_SEARCH_SPM = False
GENERATE_PLOT_CLUSTERS = True

#do not modify constants
voxel_num = VOXELS_GRID[0]*VOXELS_GRID[1]*VOXELS_GRID[2]
i = VOXELS_GRID[0]
j = VOXELS_GRID[1]
k = VOXELS_GRID[2]

#definitions
def build_voxel_coord_matrix():
	voxel_coord_matrix = np.zeros((voxel_num,3))
	cur_voxel = 0
	for ii in range (i):
		for jj in range(j):
			for kk in range(k):
				voxel_coord_matrix[cur_voxel,0] = ii
				voxel_coord_matrix[cur_voxel,1] = jj
				voxel_coord_matrix[cur_voxel,2] = kk
				cur_voxel = cur_voxel + 1

	return voxel_coord_matrix

def build_knn_matrix(voxel_coord_matrix):
	neighbors_matrix = np.zeros((K_NN-1,i*j*k))
	tree = BallTree(voxel_coord_matrix,leaf_size=2)
	for ii in range(voxel_coord_matrix.shape[0]):
		dist,ind = tree.query(voxel_coord_matrix[ii,:],k = K_NN)
		neighbors_matrix[:,ii] = ind[0,1:]
		#print str((float(ii)/float(voxel_coord_matrix.shape[0]))*100) + " %"

	return neighbors_matrix

def is_neighbor(neighbors_matrix,i,j):
	for jj in range(K_NN-1):
		if(neighbors_matrix[jj,i] == j):
			return True

	return False

def function_params(data_matrix):
	D = np.zeros((voxel_num,voxel_num))
	nnz = 0
	for ii in range(voxel_num):
		for jj in range(ii+1,voxel_num):
			D[ii,jj] = np.square(np.linalg.norm(data_matrix[ii,:] - data_matrix[jj,:]))
			if(is_neighbor(neighbors_matrix,ii,jj)):
				nnz = nnz + 1

	elem_num = ((voxel_num-1)*(voxel_num))/2
	mean_D = sum(sum(D))/elem_num
	mean_K = (2*nnz)/voxel_num
	return D,mean_D,mean_K

def J(neighbors_matrix,data_matrix):
	J_matrix = np.zeros((voxel_num,voxel_num))
	D,mean_D,mean_K = function_params(data_matrix)

	for ii in range(voxel_num):
		for jj in range(ii+1,voxel_num):
			if(is_neighbor(neighbors_matrix,ii,jj)):
				J_matrix[ii,jj] = (float(1)/float(mean_K)) * exp(float(-D[ii,jj])/float(2*mean_D*mean_D))

	return J_matrix
			
def get_m(N):
	N_max = max(N)

	m = float(Q*N_max - voxel_num)/float((Q - 1)*voxel_num)

	return m

def get_chi(m,T):
	chi = (float(voxel_num)/float(T))*float(np.var(m)) 
	return chi

def count_spins(spins):
	c = Counter(spins)

	return c.values()

def frozen_bounds(J_matrix,T,data_matrix):
	frozen_bounds_indices = {}
	for ii in range(0,J_matrix.shape[1]):
		frozen_bounds_indices[ii] = [ii]
		for jj in range(ii+1,J_matrix.shape[1]):
			if(J_matrix[ii,jj] != 0):
				if(data_matrix[ii, data_matrix.shape[1]-1] == data_matrix[jj, data_matrix.shape[1] - 1]):
					kronecker = 1
				else:
					kronecker = 0

				prob_frozen = 1 - exp((float(-1*J_matrix[ii,jj])/float(T))*kronecker)

				rnum = np.random.uniform(0, 1)
				if(rnum <= prob_frozen):
					frozen_bounds_indices[ii].append(jj)

	return frozen_bounds_indices

#stackoverflow -> macht was wir wollen
def find_clusters(aNeigh):
    def findRoot(aNode,aRoot):
        while aNode != aRoot[aNode][0]:
            aNode = aRoot[aNode][0]
        return (aNode,aRoot[aNode][1])
    myRoot = {} 
    for myNode in aNeigh.keys():
        myRoot[myNode] = (myNode,0)  
    for myI in aNeigh: 
        for myJ in aNeigh[myI]: 
            (myRoot_myI,myDepthMyI) = findRoot(myI,myRoot) 
            (myRoot_myJ,myDepthMyJ) = findRoot(myJ,myRoot) 
            if myRoot_myI != myRoot_myJ: 
                myMin = myRoot_myI
                myMax = myRoot_myJ 
                if  myDepthMyI > myDepthMyJ: 
                    myMin = myRoot_myJ
                    myMax = myRoot_myI
                myRoot[myMax] = (myMax,max(myRoot[myMin][1]+1,myRoot[myMax][1]))
                myRoot[myMin] = (myRoot[myMax][0],-1) 
    myToRet = {}
    for myI in aNeigh: 
        if myRoot[myI][0] == myI:
            myToRet[myI] = []
    for myI in aNeigh: 
        myToRet[findRoot(myI,myRoot)[0]].append(myI) 
    return myToRet  

def get_Ttrans():
 	coeff = float(4*np.log(1+np.sqrt(Q)))
 	return float(np.exp(-1/2))/float(coeff)


def in_same_cluster_cij(clusters,i,j):
	for key in clusters.keys():
		values = clusters[key]
		if((i in values) and (j in values)):
			return True

	return False

#data import
img = nib.load('diff_data.nii')
bvecs = np.loadtxt('bvecs')

print "Data import - done."

#data preprocessing
image_data = img.get_data()
n = image_data.shape[3]

bvecs_zero_ind = np.where(np.all(abs(bvecs)==0,axis=1))

voxel_coord_matrix = build_voxel_coord_matrix()
neighbors_matrix = build_knn_matrix(voxel_coord_matrix)

data_matrix = np.reshape(image_data[:i,:j,:k,:],(i*j*k,n))
data_matrix = np.delete(data_matrix,bvecs_zero_ind,axis=1)
J_matrix = J(neighbors_matrix,data_matrix)


print "Data preprocessing - done."

#algorithm

#assign random spins
spins = np.random.randint(Q, size=(data_matrix.shape[0],1))
data_matrix = np.concatenate((data_matrix,spins),axis = 1)
#SW - MC 
if(GENERATE_PLOT_SEARCH_SPM):
	chi_temp = []

	Ti = get_Ttrans() + 1 #superparamagnetic to paragmagnetic
	Tf = get_Ttrans() - float(get_Ttrans())/float(2) #ferromagnetic to superparamagnetic
	T = Ti
	temps = []
	iter_num = 0
	cluster_num = []
	while(T > Tf):
		m_steps = []
		for itr in range(M):
			frozen_bounds_indices = frozen_bounds(J_matrix,T,data_matrix)

			clusters = find_clusters(frozen_bounds_indices)
			for key in clusters.keys():
				values = clusters[key]
				data_matrix[values,data_matrix.shape[1]-1] = np.random.randint(Q)
			
			N = count_spins(data_matrix[:,data_matrix.shape[1]-1])
			m = get_m(N)
			if(itr > 9):
				m_steps.append(m)

		
		temps.append(T)
		chi = get_chi(m_steps,T)
		chi_temp.append(chi)
		cluster_num.append(len(clusters.keys()))

		T = Ti*(np.power(eta,iter_num))
		iter_num = iter_num + 1
		print iter_num

	plt.plot(temps,cluster_num,'o')
	plt.xlabel('temperature')
	plt.ylabel('cluster number')
	plt.title('Relation between cluster number and temperature')
	plt.show()

	plt.plot(temps,chi_temp,'o')
	plt.xlabel('temperature')
	plt.ylabel('chi')
	plt.title('Searching for super paramagnetic range')
	plt.show()

if(GENERATE_PLOT_CLUSTERS):
	#save C_ij in the MC, measure spin spin, get the threshold and build the clusters
	
	T_spm = 0.2 #set it to that or smth in the range (finer the range on the plot, the better we see the constant in the superparamagnetic phase)
	C_ij = np.zeros((voxel_num,voxel_num))
	spins = np.random.randint(Q, size=(data_matrix.shape[0],1))
	data_matrix = np.concatenate((data_matrix,spins),axis = 1)
	for itr in range(M):
			frozen_bounds_indices = frozen_bounds(J_matrix,T_spm,data_matrix)

			clusters = find_clusters(frozen_bounds_indices)
			for key in clusters.keys():
				values = clusters[key]
				data_matrix[values,data_matrix.shape[1]-1] = np.random.randint(Q)

			if(itr > 9):
			 	for ii in range(voxel_coord_matrix.shape[0]):
			 		for jj in range(ii+1,voxel_coord_matrix.shape[0]):
			 			if(in_same_cluster_cij(clusters,ii,jj)):
			 				C_ij[ii,jj] = C_ij[ii,jj] + 1
						


			
			print itr

	print C_ij
	print ""
	C_ij = C_ij/(M-10) # upper triangular
	print C_ij
	#sys.exit()
	C_ij = C_ij*(Q-1)
	G_ij = C_ij + np.ones((C_ij.shape[0],C_ij.shape[1]))
	G_ij = G_ij/Q # not upper triangluar- just observe the upper triangular part

	bound_clusters = dict()
	for ii in range(G_ij.shape[0]):
		bound_clusters[ii] = []
	 	for jj in range(ii+1,G_ij.shape[1]):
	 		if(G_ij[ii,jj] > 0.5):
	 			bound_clusters[ii].append(jj)
	
	new_clusters = find_clusters(bound_clusters)
	print new_clusters

	fig = plt.figure()
	ax = fig.add_subplot(111, projection='3d')
	palette = itertools.cycle(sns.color_palette())
	for key in new_clusters.keys():
		values = new_clusters[key] #all voxels in this cluster; voxel numbers
		iss = voxel_coord_matrix[values,0]
		jss = voxel_coord_matrix[values,1]
		kss = voxel_coord_matrix[values,2]

	
		c=next(palette)
		ax.scatter(iss, jss, kss,c = c, marker = ',',s = 1000)

ax.set_xlabel('i Label')
ax.set_ylabel('j Label')
ax.set_zlabel('k Label')

plt.show()
