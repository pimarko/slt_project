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
import find_clusters
import seaborn as sns

#modify constants
K_NN = 5
VOXELS_GRID = [7,7,3]
Q = 10
M = 100
eta = 0.95
GENERATE_PLOT_SEARCH_SPM = True
GENERATE_PLOT_CLUSTERS = True
coord_num = 3

#do not modify constants
voxel_num = VOXELS_GRID[0]*VOXELS_GRID[1]*VOXELS_GRID[2]
i = VOXELS_GRID[0]
j = VOXELS_GRID[1]
k = VOXELS_GRID[2]

#definitions
def build_data_matrix():	
	img = nib.load('diff_data.nii')
	bvecs = np.loadtxt('bvecs')

	image_data = img.get_data()
	n = image_data.shape[3]

	bvecs_zero_ind = np.where(np.all(abs(bvecs)==0,axis=1))

	image_data = image_data[:i,:j,:k,:]

	image_data = np.delete(image_data,bvecs_zero_ind,axis=3)

	dif_strength_num = len(bvecs) - len(bvecs_zero_ind[0])

	data_matrix = np.zeros((voxel_num,coord_num+dif_strength_num))

	spins = np.random.randint(Q, size=(voxel_num,1))

	cur_voxel = 0
	for ii in range (i):
		for jj in range(j):
			for kk in range(k):
				data_matrix[cur_voxel,0] = ii
				data_matrix[cur_voxel,1] = jj
				data_matrix[cur_voxel,2] = kk

				data_matrix[cur_voxel,3:] = image_data[ii,jj,kk,:]
				cur_voxel = cur_voxel + 1

	data_matrix = np.concatenate((data_matrix,spins),axis = 1)

	return data_matrix

def build_knn_matrix(data_matrix):
	neighbours_matrix = np.zeros((voxel_num,K_NN-1))
	tree = BallTree(data_matrix[:,0:3])
	for voxel in range(voxel_num):
		dist,ind = tree.query(data_matrix[voxel,0:3],k = K_NN)
		neighbours_matrix[voxel,:] = ind[0,1:]

	for cur_voxel in range(voxel_num):
		neighbours = neighbours_matrix[cur_voxel,:]
		for ind in range(len(neighbours)):
			neighbour = int(neighbours[ind])
			if(cur_voxel not in neighbours_matrix[neighbour,:]):
				neighbours_matrix[cur_voxel,ind] = -1

	return neighbours_matrix

def is_neighbour(neighbours_matrix,i,j):
	if(i in neighbours_matrix[j,:]):
		return True
	else:
		return False

def function_params(data_matrix,neighbours_matrix):
	D = np.zeros((voxel_num,voxel_num))
	nnz = 0
	for ii in range(voxel_num):
		for jj in range(ii+1,voxel_num):
			D[ii,jj] = np.square(np.linalg.norm(data_matrix[ii,3:-1] - data_matrix[jj,3:-1]))
			if(is_neighbour(neighbours_matrix,ii,jj)):
				nnz = nnz + 1


	elem_num = ((voxel_num-1)*(voxel_num))/2
	mean_D = sum(sum(D))/elem_num
	mean_K = (2*nnz)/voxel_num

	return D,mean_D,mean_K

def build_J_matrix(neighbours_matrix,data_matrix):
	J_matrix = np.zeros((voxel_num,voxel_num))
	D,mean_D,mean_K = function_params(data_matrix,neighbours_matrix)

	for ii in range(voxel_num):
		for jj in range(ii+1,voxel_num):
			if(is_neighbour(neighbours_matrix,ii,jj)):
				J_matrix[ii,jj] = (float(1)/float(mean_K)) * exp(float(-D[ii,jj]*D[ii,jj])/float(2*mean_D*mean_D))

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
	for ii in range(0,voxel_num):
		frozen_bounds_indices[ii] = []
		for jj in range(ii+1,voxel_num):
			if(J_matrix[ii,jj] != 0):
				if(data_matrix[ii, data_matrix.shape[1]-1] == data_matrix[jj, data_matrix.shape[1] - 1]):
					kronecker = 1
				else:
					kronecker = 0

				prob_frozen = 1 - exp((float(-1*J_matrix[ii,jj])/float(T))*kronecker)

				rnum = np.random.uniform(0, 1)
				if(rnum <= prob_frozen):
					frozen_bounds_indices[ii].append(ii)
					frozen_bounds_indices[ii].append(jj)

	return frozen_bounds_indices

def get_Ttrans():
 	coeff = float(4*np.log(1+np.sqrt(Q)))
 	return float(np.exp(float(-1/2)))/float(coeff)


def in_same_cluster_cij(clusters,i,j):
	for key in clusters.keys():
		values = clusters[key]
		if((i in values) and (j in values)):
			return True

	return False


#data preprocessing
data_matrix = build_data_matrix()
neighbours_matrix = build_knn_matrix(data_matrix)
J_matrix = build_J_matrix(neighbours_matrix,data_matrix)

print "Data preprocessing - done."

#algorithm
if(GENERATE_PLOT_SEARCH_SPM):
	chi_temp = []

	Ti = get_Ttrans()*10
	Tf = get_Ttrans()*0.1
	T = Ti
	temps = []
	iter_num = 0
	cluster_num = []
	m_means = []
	while(T > Tf):
		m_steps = []
		for itr in range(M):
			frozen_bounds_indices = frozen_bounds(J_matrix,T,data_matrix)

			clusters = find_clusters.find_clusters(frozen_bounds_indices)
			for key in clusters.keys():
				values = clusters[key]
				data_matrix[values,-1] = np.random.randint(Q)
			
			N = count_spins(data_matrix[:,-1])
			m = get_m(N)
			if(itr > 9):
				m_steps.append(m)

		
		temps.append(T)
		chi = get_chi(m_steps,T)
		m_means.append(np.mean(m_steps))
		chi_temp.append(chi)
		cluster_num.append(len(clusters.keys()))

		T = Ti*(np.power(eta,iter_num))
		iter_num = iter_num + 1
		print iter_num

	plt.plot(temps,cluster_num,'-b')
	plt.xlabel('temperature')
	plt.ylabel('cluster number')
	plt.title('Relation between cluster number and temperature')
	plt.show()

	plt.plot(temps,chi_temp,'-b')
	plt.xlabel('temperature')
	plt.ylabel('chi')
	plt.title('Searching for super paramagnetic range')
	plt.show()

	plt.plot(temps,m_means,'-b')
	plt.xlabel('temperature')
	plt.ylabel('<m>')
	plt.title('Relation between <m> and temperature')
	plt.show()

if(GENERATE_PLOT_CLUSTERS):
	index_temp = chi_temp.index(max(chi_temp)) + 1
	T_spm = temps[index_temp]
	print T_spm
	C_ij = np.zeros((voxel_num,voxel_num))
	spins = np.random.randint(Q, size=(voxel_num,1))
	spins = np.reshape(spins,(voxel_num,))
	data_matrix[:,-1] = spins
	for itr in range(M):
			frozen_bounds_indices = frozen_bounds(J_matrix,T_spm,data_matrix)

			clusters = find_clusters.find_clusters(frozen_bounds_indices)
			for key in clusters.keys():
				values = clusters[key]
				data_matrix[values,-1] = np.random.randint(Q)

			if(itr > 9):
			 	for ii in range(voxel_num):
			 		for jj in range(ii+1,voxel_num):
			 			if(in_same_cluster_cij(clusters,ii,jj)):
			 				C_ij[ii,jj] = C_ij[ii,jj] + 1
						


			
			print itr

	C_ij = C_ij/(M-10) 
	C_ij = C_ij*(Q-1)
	G_ij = C_ij + np.ones((C_ij.shape[0],C_ij.shape[1]))
	G_ij = G_ij/Q

	plt.hist(np.reshape(G_ij,(G_ij.shape[0]*G_ij.shape[1],1)))
	plt.xlabel('G_ij')
	plt.ylabel('#')
	plt.title('temperature ' + str(T_spm))
	plt.show()

	bound_clusters = dict()
	for ii in range(G_ij.shape[0]):
		bound_clusters[ii] = []
	 	for jj in range(ii+1,G_ij.shape[1]):
	 		if(G_ij[ii,jj] > 0.5):
	 			bound_clusters[ii].append(ii)
	 			bound_clusters[ii].append(jj)
	
	new_clusters = find_clusters.find_clusters(bound_clusters)

	fig = plt.figure()
	ax = fig.add_subplot(111, projection='3d')
	palette = itertools.cycle(sns.color_palette())
	for key in new_clusters.keys():
		values = new_clusters[key]
		iss = data_matrix[values,0]
		jss = data_matrix[values,1]
		kss = data_matrix[values,2]

	
		c = next(palette)
		ax.scatter(iss, jss, kss,c = c, marker = ',',s = 1000)

	ax.set_xlabel('i Label')
	ax.set_ylabel('j Label')
	ax.set_zlabel('k Label')

	plt.show()
