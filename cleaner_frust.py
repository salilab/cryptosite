import sys

data = open(sys.argv[-1])
D = data.readlines()
data.close()

out = open(sys.argv[-1],'w')
for d in D:
	if d[:4]=='ATOM' or d[:3]=='TER': out.write(d[:21]+'A'+d[22:])
out.close()



