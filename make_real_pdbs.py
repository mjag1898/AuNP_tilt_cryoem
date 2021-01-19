from scipy.io import loadmat
import numpy as np
from tqdm import tqdm

structs_file = 'real_structs.mat'
coords = loadmat(structs_file)['structs'] 

gold_lines = open('golds.pdb', 'r').read().splitlines()
gold_coords = np.array([line[31:55].split() for line in gold_lines]).astype(np.float32)
gold_coords -= gold_coords.mean(axis=0)

full_coords = np.zeros((coords.shape[2], coords.shape[0], gold_coords.shape[0], 3))
for i in range(full_coords.shape[0]):
	for j in range(full_coords.shape[1]):
		full_coords[i, j, :, :] = gold_coords + coords[j, :, i]

chains = ['A', 'B', 'C', 'D', 'E', 'F', 'G']
for i in tqdm(range(full_coords.shape[0])):
	new_lines = []
	for j in range(full_coords.shape[1]):
		new_gold_lines = gold_lines.copy()
		for k in range(full_coords.shape[2]):
			new_gold_lines[k] = (
				new_gold_lines[k][:21] +
				chains[j] +
				new_gold_lines[k][22:30] + 
				"{:8.3f}".format(full_coords[i, j, k, 0]) + 
				"{:8.3f}".format(full_coords[i, j, k, 1]) +
				"{:8.3f}".format(full_coords[i, j, k, 2]) +
				new_gold_lines[k][54:]
			)
		new_lines += new_gold_lines
	open('real_struct_'+str(i+1)+'.pdb', 'w').write("\n".join(new_lines))

for i in tqdm(range(full_coords.shape[0])):
	new_lines = []
	for j in range(full_coords.shape[1]):
		gold_line = gold_lines[0]
		com_coord = full_coords[i, j].mean(axis=0)
		gold_line = (
			gold_line[:21] +
			chains[j] +
			gold_line[22:30] +
			"{:8.3f}".format(com_coord[0]) + 
			"{:8.3f}".format(com_coord[1]) +
			"{:8.3f}".format(com_coord[2]) +
			gold_line[54:]
		)
		new_lines.append(gold_line)
	open('real_struct_com_'+str(i+1)+'.pdb', 'w').write("\n".join(new_lines))

