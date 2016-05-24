print(__doc__)

from time import time

# these two lines have to be inserted here
#import matplotlib
#matplotlib.use('Agg')

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
K_NN = 26
VOXELS_GRID = [5,5,5]
Q = 10
M = 20
eta = 0.95
GENERATE_PLOT_SEARCH_SPM = True
GENERATE_PLOT_CLUSTERS = True
coord_num = 3
thres = 10

#do not modify constants
voxel_num = VOXELS_GRID[0]*VOXELS_GRID[1]*VOXELS_GRID[2]
init_i = 135
init_j = 76
init_k = 74
i = VOXELS_GRID[0] + init_i
j = VOXELS_GRID[1] + init_j
k = VOXELS_GRID[2] + init_k

i_step = VOXELS_GRID[0] 
j_step = VOXELS_GRID[1] 
k_step = VOXELS_GRID[2] 


#definitions
def build_data_matrix():	
	img = nib.load('diff_data.nii')
	bvecs = np.loadtxt('bvecs')

	image_data = img.get_data()
	n = image_data.shape[3]

	bvecs_zero_ind = np.where(np.all(abs(bvecs)==0,axis=1))

	image_data = image_data[init_i:i,init_j:j,init_k:k,:]

	image_data = np.delete(image_data,bvecs_zero_ind,axis=3)

	dif_strength_num = len(bvecs) - len(bvecs_zero_ind[0])

	data_matrix = np.zeros((voxel_num,coord_num+dif_strength_num))

	spins = np.random.randint(Q, size=(voxel_num,1))

	cur_voxel = 0
	for ii in range (i_step):
		for jj in range(j_step):
			for kk in range(k_step):
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
	D_normal = np.zeros((voxel_num,voxel_num))
	nnz = 0
	for ii in range(voxel_num):
		for jj in range(ii+1,voxel_num):
			D[ii,jj] = np.square(np.linalg.norm(data_matrix[ii,3:data_matrix.shape[1]-2] - data_matrix[jj,3:data_matrix.shape[1]-2]))
			D_normal[ii,jj] = np.linalg.norm(data_matrix[ii,3:data_matrix.shape[1]-2] - data_matrix[jj,3:data_matrix.shape[1]-2])
			if(is_neighbour(neighbours_matrix,ii,jj)):
				nnz = nnz + 1


	elem_num = float((voxel_num-1)*(voxel_num))/float(2)
	mean_D = float(sum(sum(D)))/float(elem_num)
	mean_K = float(2*nnz)/float(voxel_num)
	mean_D_normal = float(sum(sum(D_normal)))/float(elem_num)

	return D,mean_D,mean_K,mean_D_normal

def build_J_matrix(neighbours_matrix,data_matrix):
	J_matrix = np.zeros((voxel_num,voxel_num))
	D,mean_D,mean_K,mean_D_normal = function_params(data_matrix,neighbours_matrix)

	for ii in range(voxel_num):
		for jj in range(ii+1,voxel_num):
			if(is_neighbour(neighbours_matrix,ii,jj)):
				J_matrix[ii,jj] = (float(1)/float(mean_K)) * exp(float(-D[ii,jj])/float(2*mean_D_normal*mean_D_normal))

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
		frozen_bounds_indices[ii] = [ii]
		for jj in range(ii+1,voxel_num):
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

def get_Ttrans(data_matrix,neighbours_matrix):
	D,mean_D,mean_K,mean_D_normal = function_params(data_matrix,neighbours_matrix)
 	coeff = float(4*float(np.log10(1+float(np.sqrt(Q)))))
 	exponent = -1*float(mean_D)/float(2*mean_D_normal*mean_D_normal)
 	return float(np.exp(float(exponent)))/float(coeff)

def in_same_cluster_cij(clusters,ii,jj):
	for key in clusters.keys():
		values = clusters[key]
		if((ii in values) and (jj in values)):
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
	ttrans = get_Ttrans(data_matrix,neighbours_matrix)
	Ti = ttrans*4
	Tf = float(ttrans)*0.1
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
	#plt.savefig('plot_cluster_num_temps.png')
	plt.show()

	index_temp = chi_temp.index(max(chi_temp)) - 1
	T_spm = temps[index_temp]
	print T_spm

	plt.plot(temps,chi_temp,'-b',temps[index_temp],chi_temp[index_temp],'-o')
	plt.xlabel('temperature')
	plt.ylabel('chi')
	plt.title('Searching for super paramagnetic range')
	#plt.savefig('plot_chi_temps.png')
	plt.show()
	
	plt.plot(temps,m_means,'-b')
	plt.xlabel('temperature')
	plt.ylabel('<m>')
	plt.title('Relation between <m> and temperature')
	#plt.savefig('plot_m_temps.png')
	plt.show()
	
if(GENERATE_PLOT_CLUSTERS):
	C_ij = np.zeros((voxel_num,voxel_num))
	G_ij = np.zeros((voxel_num,voxel_num))
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

	for ii in range(voxel_num):
		for jj in range(ii+1,voxel_num):
			C_ij[ii,jj] = float(C_ij[ii,jj])/float(M-10) 
			G_ij[ii,jj] = float((C_ij[ii,jj]*(Q-1)) + 1)/float(Q)

	array_Gij = []
	for ii in range(G_ij.shape[0]):
		for jj in range(ii+1,G_ij.shape[1]):
			array_Gij.append(G_ij[ii,jj])

	plt.hist(array_Gij)
	plt.xlabel('G_ij')
	plt.ylabel('#')
	plt.title('temperature ' + str(T_spm))
	#plt.savefig('plot_gij_hist.png')
	plt.show()
	
	bound_clusters = dict()
	for ii in range(G_ij.shape[0]):
		bound_clusters[ii] = [ii]
	 	for jj in range(ii+1,G_ij.shape[1]):
	 		if(G_ij[ii,jj] > 0.5):
	 			bound_clusters[ii].append(jj)
	
	new_clusters = find_clusters.find_clusters(bound_clusters)

	print new_clusters
	print "keys:" + str(new_clusters.keys())

	fig = plt.figure()
	ax = fig.add_subplot(111, projection='3d')
	mmax = np.max(new_clusters.keys())
	for key in new_clusters.keys():
		values = new_clusters[key]
		if(len(values) > thres):
			iss = data_matrix[values,0]
			jss = data_matrix[values,1]
			kss = data_matrix[values,2]

			ax.scatter(iss, jss, kss,c = [[key/float(mmax),key/float(mmax),key/float(mmax)]], marker = ',',s = 2000)


	ax.set_xlabel('i')
	ax.set_ylabel('j')
	ax.set_zlabel('k')

	#plt.savefig('plot_final_clusters.png')
	plt.show()