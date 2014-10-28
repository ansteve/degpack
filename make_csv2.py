import os

def make_csv(outdir, filter_names, filename, clsList, ntop):
	if(outdir != None):
		if not os.path.exists(outdir):
			os.makedirs(outdir)

	filter_names = [tup[0] for tup in filter_names]
	filter_names = filter_names[0:ntop]

	fin = open(filename, 'r')
	outfile_name = 'out'
	
	fout = open(outdir+'/'+outfile_name+'_top'+str(ntop)+'.csv', 'w')
	fcls = fin.readline().strip().split()
	
	outDict = {}
	for line in fin:
		spl = line.strip().split()
		outDict[spl[0]] = spl[1:]
	fin.close()

	for gene in filter_names:
		fout.write(gene+',')
	fout.write('class\n')
	
	for i in range(len(clsList)):
		for gene in filter_names:
			fout.write(outDict[gene][i]+',')
		fout.write('c'+clsList[i])
		fout.write('\n')
	fout.close()

