#! python2.7
##################  Cell Type and Gene Expression Simulation  ##################
import numpy as np
import random
from bioinf import express
from hamming import print_hamming
random.seed()

##################  Select Parameters  ##################
cells = 100 ;	genes = 10000	;  cell_types = 5  ;	gene_sets = 25  ;    mode = 'uniform'

cell_profile = np.array([ \
['o','o','o','h','h','h','h','l','l','l','o','o','o','o','o','o','o','o','o','o','o','o','o','o','o'],\
['o','o','o','h','h','h','h','l','l','l','l','l','l','o','o','o','o','o','o','o','o','o','o','o','o'],\
['o','o','o','o','o','o','o','o','o','o','o','o','o','o','o','o','o','o','o','o','o','o','o','o','o'],\
['o','o','o','o','o','o','o','o','o','o','o','o','o','o','o','o','o','o','o','o','o','o','o','o','o'],\
['o','o','o','o','o','o','o','o','o','o','o','o','o','o','o','o','o','o','o','o','o','o','o','o','o']])


cell_type_distribution = [20, 20, 20, 20, 20];
	
gene_set_distribution = [400, 400, 400, 400, 400,\
			 400, 400, 400, 400, 400,\
			 400, 400, 400, 400, 400,\
		  	 400, 400, 400, 400, 400,\
			 400, 400, 400, 400, 400]

cell_types_to_sample = np.array([0,1,2,3])

#Randomly assign based on cell type expression level values to each gene.
a = np.zeros((cells, genes),dtype='float')
I,J = a.shape
cell_profile_expanded = np.repeat(np.repeat(cell_profile,cell_type_distribution,axis=0),gene_set_distribution,axis=1)
for i in range(I):
	for j in range(J):#consider xrange()
		a[i][j] = express(mode,cell_profile_expanded[i][j])
		

#Sample 4 cell types with replacement from available cell types. Randomly choose a cell(s) from each cell type.

#cell_types_to_sample = np.array(random.sample(4*range(cell_types),4))
types_, frequency = np.unique(cell_types_to_sample, return_counts=True); cells_sampled=[];
cumsum = np.cumsum(cell_type_distribution)
for i in range(len(types_)):
	type_, boundary= cell_type_distribution[types_[i]], cumsum[types_[i]]
	cells_sampled += list(random.sample(range(type_),frequency[i])+boundary-type_)	

b = np.zeros((3,4,genes),dtype='float'); permutations = [[0,1,2,3],[0,3,1,2],[0,1,3,2]]
b[0,:,:] = a[cells_sampled]
b[1,:,:] = b[0,:,:][permutations[1]]
b[2,:,:] = b[0,:,:][permutations[2]]


#Compute 3 permutations of Hamming distance matrix for each pair of cells
c = np.zeros((3,4,4),dtype='float')
K,I,J = c.shape
for k in range(K):
	for i in range(I):
		for j in range(J):
			c[k][i][j] = sum(abs(b[k][i][:]-b[k][j][:]))



#Compute Hamming distance box plots for each of the permutations. Discard permutations that result in negative variables.
print('Cells:%i\nGenes:%i\nExpression:random_%s_value_assigned_to_each_gene\nHamming matrix:'%(cells,genes,mode))
for k in range(K):
	H1,H2,H3=c[k][0][1:]
	H4,H5=c[k][1][2:]
	H6=c[k][2][3]
	A=int(0.5*(H1+H2-H4))
	B=int(0.5*(H1+H5-H3))
	C=int(0.5*(H2+H6-H3))
	D=int(0.5*(H5+H6-H4))
	E=int(0.5*(H3+H4-H2-H5))
	F=int(0.5*(H3+H4-H1-H6))
	coordinates = np.array([A,B,C,D,E,F])/(max(A,B,C,D,E,F)/10); delta=min(E,F)/float(max(E,F))
	X,Y,U,V = cell_types_to_sample[permutations[k]]
	if (A<0 or B<0 or C<0 or D<0 or E<0 or F<0):
		print('\tPermutation %i: negative coordinates'%k)
		print('\tcoordinates= [%i,%i,%i,%i,%i,%i]'%(A,B,C,D,E,F))
		print('\tcell_type_at_vertices= [%i,%i,%i,%i]\n'%(X,Y,U,V))
		continue
	else:
		print('\tPermutation %i:'%k)
		print('\tdelta= %f'%delta)
		print('\tcoordinates= [%i,%i,%i,%i,%i,%i]'%(A,B,C,D,E,F))
		print('\tcell_type_at_vertices= [%i,%i,%i,%i]\n'%(X,Y,U,V))
		print_hamming(coordinates)



#Compute completely random Hamming distance matrix
#d = np.zeros((4,4),dtype='int')
#I,J = d.shape
#for i in range(I):
#	for j in range(J):
#		if j>i:
#			d[i][j] = d[j][i] = random.randint(0,10000)

#Sample 4 cells without replacement from total cell population. Determine cell types from sample cell index.
#cells_sampled = random.sample(range(cells),4);  cell_types_sampled= []
#for cell_sampled in cells_sampled:
#	j = 0; cumsum = 0
#	while j < cell_types:
#		cumsum += cell_type_distribution[j]
#		if cell_sampled < cumsum:
#			cell_types_sampled.append(j)
#			break
#		j+=1 













