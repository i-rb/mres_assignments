import pandas as pd
import math as mt
import numpy as np
import networkx as nx
import copy
import scipy.sparse.linalg
import glee

# outputs path
embedding_network_path = "/Users/ivanrendobarreiro/Documents/academico/mres/s2/ML/matrix_to_impute_NET.csv"
real_data_formatted_path = "/Users/ivanrendobarreiro/Documents/academico/mres/s2/ML/matrix_real.csv"

# input path
path_to_original_data = "/Users/ivanrendobarreiro/Documents/academico/mres/s2/ML/ia-workplace-contacts/ia-workplace-contacts.txt"
data = pd.read_csv(path_to_original_data,delimiter=",",header=0)

# transform time to days
df["t"] = np.floor(df["t"] / (6*86400))

# create nodes by groups
nodes = np.unique([df["j"], df["i"]])
nodes_norm = np.array(list(range(len(nodes))))

df2 = df.copy()
node_mapping = {node: norm for node, norm in zip(nodes, nodes_norm)}
df2["i"] = df["i"].map(node_mapping)
df2["j"] = df["j"].map(node_mapping)
df2["index"] = range(len(df2))
df2 = df2.set_index("index")

conditions = [(df2['i'].between(0, 10)) & (df2['j'].between(0, 10)),
              (df2['i'].between(11, 20)) & (df2['j'].between(11, 20)),
              (df2['i'].between(21, 30)) & (df2['j'].between(21, 30)),
              (df2['i'].between(31, 40)) & (df2['j'].between(31, 40)),
              (df2['i'].between(41, 50)) & (df2['j'].between(41, 50)),
              (df2['i'].between(51, 60)) & (df2['j'].between(51, 60)),
              (df2['i'].between(61, 70)) & (df2['j'].between(61, 70)),
              (df2['i'].between(71, 80)) & (df2['j'].between(71, 80)),
              (df2['i'].between(81, 90)) & (df2['j'].between(81, 90))]

df2['label'] = np.select(conditions, ['0','1', '2', '3','4','5','6','7','8'], default=None)

filtered_df = df2.dropna(subset=['label'])

# from here, complete the 9x2 matrix with the adj matrices of each group at each t

mat = []
for lab in range(9):
    mat_lab = []
    for t in range(2):
        tempdf = filtered_df[(filtered_df["label"]==str(lab))&(filtered_df["t"]==t)]
        adj_matrix = np.zeros((10, 10))
        for _, row in tempdf.iterrows():
            i, j = row['i'], row['j']
            adj_matrix[i%10, j%10] = 1
            adj_matrix[j%10, i%10] = 1
        mat_lab.append(adj_matrix)
    mat.append(mat_lab)
    
mat_real = copy.deepcopy(mat)


# create nans

mat[4][1] = np.nan
mat[3][0] = np.nan
mat[6][1] = np.nan
mat[7][1] = np.nan


init2 = np.zeros((1,20))
for lab in range(9):
    if np.isnan(mat[lab][0]).any():
        mat[lab][0] = np.full((10, 10), np.nan)
    if np.isnan(mat[lab][1]).any():
        mat[lab][1] = np.full((10, 10), np.nan)
    hor2 = np.hstack((mat_real[lab][0], mat_real[lab][1]))
    init2 = np.vstack((init2, hor2))
init2 = np.delete(init2, 0, axis=0)


# using GLEE to embed the matrix

init = np.zeros((1,2))
for i in mat:
    G_temp = nx.from_numpy_matrix(i[0])
    E0 = eigenmaps(G_temp, dim=1, method='glee', return_vals=False) # this is the GLEE command
    G_temp = nx.from_numpy_matrix(i[1])   
    E1 = eigenmaps(G_temp, dim=1, method='glee', return_vals=False) # this is the GLEE command
    hor = np.hstack((E0, E1))
    init = np.vstack((init, hor))
init = np.delete(init, 0, axis=0)

np.savetxt(embedding_network_path, init, delimiter=';')
np.savetxt(real_data_formatted_path, init2, delimiter=';')