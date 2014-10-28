import random
import math
import numpy as np
import matplotlib.pyplot as plt
pseudo_cnt = 0.0000001

def mi(vec1, vec2):
	len1 = len(vec1)
	len2 = len(vec2)
	#p1 = {0: 0.5, 1: 0.5}
	p1 = {}
	p2 = {}
	for elem in vec1:
		if not elem in p1:
			p1[elem] = 0.0
		p1[elem] += 1.0
	for key in p1.keys():
		p1[key] = (p1[key] / len1)
	for elem in vec2:
		if not elem in p2:
			p2[elem] = 0.0
		p2[elem] += 1.0
	for key in p2.keys():
		p2[key] = (p2[key] / len2)
	
	p12 = {}
	for k1 in p1.keys():
		for k2 in p2.keys():
			k12 = str(k1)+','+str(k2)
			if not k12 in p12:
				p12[k12] = pseudo_cnt
	for i in range(len1):
		key = str(vec1[i])+','+str(vec2[i])
		p12[key] += 1.0
	#print len1
	#print p12
	for key in p12.keys():
		p12[key] = p12[key] / len1
	
	#print p12
	mi = 0.0
	for k1 in p1.keys():
		for k2 in p2.keys():
			k12 = str(k1)+','+str(k2)
			mi += p12[k12]*math.log((p12[k12]/(p1[k1]*p2[k2])),2)
	#print mi
	return mi


def rand_permut(clsList, bin_num):
	numSample = len(clsList)
	cls2idx = {}
	idx=0
	for cls in clsList:
		if not cls in cls2idx:
			cls2idx[cls]=idx
			idx+=1
	
	ranks = []
	#for cls in clsList:
	#	ranks.append(cls2idx[cls])
	for i in range(len(clsList)):
		r = (i / bin_num) + 1
		ranks.append(r)

	res = []
	for i in range(100000):
		res.append(mi(clsList, ranks))
		random.shuffle(ranks)
		#print ranks
		#print res[-1]
	hist, bin_edges = np.histogram(res)
	print hist, bin_edges
	plt.hist(res)
	plt.savefig('temp.png')
	return hist, bin_edges

rand_permut([1,1,1,1,1,2,2,2,2,2,3,3,3,3,3],3)
