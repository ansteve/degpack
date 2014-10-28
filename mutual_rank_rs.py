#!/usr/bin/env python
import math
import sys
from my_argparse import *
from make_csv2 import *
import random
import numpy as np

pseudo_cnt = 0.000000001

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

def discretize_values(values, bin_num):
	sorted_vals = sorted(values)
	length = len(sorted_vals)
	bin_edges = []
	for i in range(bin_num):
		bin_edges.append(sorted_vals[length * i / bin_num])
	bin_edges.append(float('inf'))
	#print bin_edges
	ranks = []
	for val in values:
		for i in range(bin_num):
			if bin_edges[i] <= val and val < bin_edges[i+1]:
				ranks.append(i)
				break
	#print ranks
	return ranks

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
	hist, bin_edges = np.histogram(res, bins=50)
	return hist, bin_edges

def get_pvalue(hist, bin_edges, score):
	sum_hist = sum(hist)
	pre_pval = 0
	for i in range(0, len(hist)):
		if bin_edges[i+1] >= score:
			pre_pval += hist[i]
	pvalue = float(pre_pval) / sum_hist
	return pvalue

# run_mutual : main function. calculate the mutual information of every gene.
# return form : each line has (gene_id, score, p_value) 
def run_mutual(input_matrix, geneList, clsList, bin_num):
	ret_lines = []
	hist, bin_edges = rand_permut(clsList, bin_num)
	for li in range(len(input_matrix)):
		values = input_matrix[li]
		id = geneList[li]
		disc_values = discretize_values(values, bin_num)
		score = mi(clsList, disc_values)
		pvalue = get_pvalue(hist, bin_edges, score)
		ret_lines += [(id, score, pvalue)]
		#ret_lines += [(id, score, 0)]

	return ret_lines

# check_num_samples : if the numbers of samples are same, return true
#				else, return flase
# return form : true or false, minimum number of samples
def check_num_samples(clsList):
	clsHash = {} # class to the number of samples
	for cls in clsList:
		if not cls in clsHash:
			clsHash[cls] = 0
		clsHash[cls] += 1
	minNum = 100000000
	maxNum = 0
	for cls in clsHash.keys():
		if clsHash[cls] <= minNum:
			minNum = clsHash[cls]
		if clsHash[cls] >= maxNum:
			maxNum = clsHash[cls]
	if minNum == maxNum:
		return True, minNum
	else:	
		return False, minNum

# random_sampling : make new fmat to make the numbers of samples 
# return form : fmat class
def random_sampling(oldFmat, numSample):
	clsHash = {} # class to index list
	totalClsHash = {} # index to class
	newClsHash = {} # class to index list
	for i in range(len(oldFmat.clsList)):
		cls = oldFmat.clsList[i]
		if not cls in clsHash:
			clsHash[cls] = []
		clsHash[cls].append(i)
	for cls in clsHash.keys():
		newClsHash[cls] = random.sample(clsHash[cls], numSample)
		for i in newClsHash[cls]:
			totalClsHash[i] = cls
	indexList = totalClsHash.keys()
	indexList.sort()
	# set newFmat
	newFmat = fin_matrix()
	newFmat.geneList = oldFmat.geneList
	newFmat.clsList = []
	for i in indexList:
		newFmat.clsList.append(totalClsHash[i])
	newFmat.input_matrix = [[] for j in range(len(oldFmat.input_matrix))]
	for j in range(len(oldFmat.input_matrix)):
		for i in indexList:
			newFmat.input_matrix[j].append(oldFmat.input_matrix[j][i])
	newFmat.nSample = len(indexList)
	newFmat.nClass = len(clsHash.keys())
	return newFmat

if __name__=='__main__':
	info = arg_parsing()
	fmat = fin_parsing(info.fin_name, info.cond)
	#binNum = int(math.log(len(fmat.clsList),2)+1)
	binNum = fmat.nClass
	
	# if the numbers of samples in other groups are not same, do random sampling to make the numbers same.
	isSame, minNumSample = check_num_samples(fmat.clsList)
	if(isSame==False):
		totalGeneScores = {} # store gene scores for each iteration
		# initialize totalGeneScores.
		for gene in fmat.geneList:
			totalGeneScores[gene] = [] 
		for i in range(2): # repeat 10 times.
			newFmat = random_sampling(fmat, minNumSample)
			output_tuples = run_mutual(newFmat.input_matrix, newFmat.geneList, newFmat.clsList, binNum)
			for tup in output_tuples:
				(gene, score) = (tup[0], tup[1])
				totalGeneScores[gene].append(score)
		# calculate a new score by averaging the old scores.
		result_tuples = []
		for gene in totalGeneScores.keys():
			totalGeneScores[gene] = np.mean(totalGeneScores[gene])
			result_tuples.append((gene,totalGeneScores[gene]))
		result_tuples.sort(key=lambda tuple: tuple[1], reverse=True)
		# print the results
		i=1
		for tup in result_tuples[0:info.ntop]:
			print str(i)+'\t'+tup[0]+'\t'+str(tup[1])
			i+=1
	else:
		output_tuples = run_mutual(fmat.input_matrix, fmat.geneList, fmat.clsList, binNum)
		output_tuples.sort(key=lambda tuple: tuple[1], reverse=True)
		#print 'Order\tName\tScore\tP-value'
		i=1
		for tup in output_tuples[0:info.ntop]:
			print str(i)+'\t'+tup[0]+'\t'+str(tup[1])+'\t'+str(tup[2])
			#print tup[0]
			i+=1

	#if(info.outdir != None):
	#	make_csv(info.outdir, output_tuples, info.fin_name, fmat.clsList, 50)
