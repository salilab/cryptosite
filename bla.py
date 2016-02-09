import glob

files = glob.glob('tt*/output/*/pred_dE0.1rAS1000/*/allosmod.py')
print len(files)
P = {}
for fil in files:
	pdb = fil.split('/')[2]
	cp = fil.split('/')[4]	
	if pdb not in P: P[pdb] = [cp]
	else: P[pdb].append(cp)

import pickle
out = open('check.pkl','w')
pickle.dump(P,out)
out.close()

