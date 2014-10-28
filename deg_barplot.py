import matplotlib.pyplot as plt
import sys

fin = open(sys.argv[1], 'r')
deg = sys.argv[2]

fin.readline() # first line is header
for line in fin:
	spl = line.strip().split()
	if(deg == spl[0]):
		pspl = [float(x) for x in spl[1:]]
		break

ax = plt.subplot(111)
ax.bar([1,2,3,4,5,6,7,8,9], pspl, width=0.1)

plt.show()
