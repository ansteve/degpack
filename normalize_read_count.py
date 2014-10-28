#!/usr/bin/env python
import sys
import numpy as np

## argument processing
if(len(sys.argv)!=4):
	print "Usage: ./this_code.py input_table.txt [depth | median] output.txt"
	sys.exit(0)
if(sys.argv[2]!='depth' and sys.argv[2]!='median'):
	print "Error: Second argument should be rpkm or median."
	sys.exit(0)
method = sys.argv[2]

## When the normalization method uses sequencing depth,
## the factor is a library size.
fin = open(sys.argv[1], 'r')
header = fin.readline().split()
factor = []
if(method=='depth'):
	factor = [0 for i in range(1, len(header))]
	for line in fin:
		spl = line.strip().split()
		for i in range(1, len(spl)):
			factor[i-1] += float(spl[i])
	fin.close()

## When the normalization method uses median read count value,
## the factor is a median value.
elif(method=='median'):
	factor = [[] for i in range(1, len(header))]
	for line in fin:
		spl = line.strip().split()
		for i in range(1, len(spl)):
			factor[i-1].append(float(spl[i]))
	for i in range(len(factor)):
		factor[i] = np.median(factor[i])
	fin.close()

## post processing factor
print factor
max_factor = max(factor)
factor = [max_factor/x for x in factor]
print factor

## print normalized expression values
fin = open(sys.argv[1], 'r')
fout = open(sys.argv[3], 'w')
spl=fin.readline().strip().split()
fout.write('\t'.join(spl)+'\n')
for line in fin:
	spl=line.strip().split()
	fout.write(spl[0])
	for i in range(1, len(spl)):
		spl[i] = str(round(float(spl[i])*factor[i-1],4))
		fout.write('\t'+spl[i])
	fout.write('\n')

