#!/usr/bin/env python
import sys

ntop = 50

poiseq = open(sys.argv[1],'r')
fout = open(sys.argv[1]+'.50.csv', 'w')
#fin = open('simulation/simul.txt', 'r')
fin = open(sys.argv[2], 'r')
fcond = open(sys.argv[3], 'r')

glist = []
poiseq.readline()
for line in poiseq:
	gname = line.split()[1].strip('"')
	glist.append(gname)
	if(len(glist)==ntop):
		break

print glist

outmat = []
inmat = {}
fin.readline() # first line is header
for line in fin:
	spl = line.strip().split()
	inmat[spl[0]] = spl

for id in glist:
	outmat.append(inmat[id])

cline = ['class']
line = fcond.readline()
for cls in line.strip().split():
	cline.append('c'+cls)

outmat.append(cline)
print outmat
for i in range(len(outmat[0])):
	fout.write(outmat[0][i])
	for j in range(1, len(outmat)):
		fout.write(','+outmat[j][i])
	fout.write('\n')

