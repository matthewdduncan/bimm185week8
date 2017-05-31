##note that certain input files were modified very slightly to prevent major outliers generated from bad data

import matplotlib.pyplot as plt
import numpy
from scipy.stats import gaussian_kde

##get only strong and confirmed, split by strand
with open('OperonSet.txt') as oset:
	with open('operonf.txt', 'w') as outf:
		for line in oset:
			if '#' in line[0:1]:
				continue
			if 'Confirmed' in line or 'Strong' in line:
				outf.write(line)

## create a list of all operons				
with open('operonf.txt') as forward:
	operonsf = []
	for line in forward:
		k = line.split('\t')
		operonsf.append(k[0:6])
	operonsf.sort(key=lambda x:int(x[1]))

##create a dict with allgenes
with open('GeneProductSet.txt') as genes:
	genetable = {}
	for line in genes:
		if '#' in line[0:1] or 'Phantom Gene' in line:
			continue
		k = line.split('\t')
		if k[3] != '':
			genetable[k[1]] = [int(k[3]), int(k[4])]

##simple function for calculating the distance between genes			
def distance(a, b, genetable):
	if genetable[a][0] > genetable[b][0]:
		return genetable[a][0] - genetable[b][1] + 1
	else:
		return genetable[b][0] - genetable[a][1] + 1

slgenes = []

#sorted list of genes
for gene in genetable.keys():
	slgenes.append([gene, genetable[gene][0], genetable[gene][1]])
	
slgenes.sort(key=lambda x:x[1])

##distance list for genes in operon		
opdist = []
ndist = []

#populate the distance list in operon
for op in operonsf:
	prevg = ''
	gio = op[5].split(',')
	for g in gio:
		if prevg != '':
			opdist.append(distance(prevg, g, genetable))
		
		for pos,gene in enumerate(slgenes):
			if g in gene:
				if slgenes[pos-1][0] not in gio:
					if distance(g, slgenes[pos-1][0], genetable) < 50000:
						ndist.append(distance(g, slgenes[pos-1][0], genetable))
				if slgenes[pos+1][0] not in gio:
					if distance(g, slgenes[pos+1][0], genetable) < 50000:
						ndist.append(distance(g, slgenes[pos+1][0], genetable))
		prevg = g
	
		
#sort distance lists		
opdist.sort()
ndist.sort()

#Plot density functions for both hypotheses
axes = plt.gca()
axes.set_xlim([-400,1000])

density = gaussian_kde(ndist)
xs = numpy.linspace(-1000,1500,1000)
density.covariance_factor = lambda : .8 #.25 on report
density._compute_covariance()
plt.plot(xs,density(xs))

density2 = gaussian_kde(opdist)
density2.covariance_factor = lambda : .8 #.25 on report
density2._compute_covariance()
plt.plot(xs,density2(xs))

plt.legend(['Negative Control', 'Postive Control'])

plt.xlabel('Gene Distance')
plt.ylabel('Density')

plt.savefig('KDE.png')

plt.close('all')

#evaluates the likelyhood of h1, given both density functions, the distance, and the p(h1)
def evaluateit(d1, d2, dist, factor):
	return (d2.evaluate(dist)*factor)/(d2.evaluate(dist)*factor + d1.evaluate(dist)*(1-factor))

fr = .6
info = []
info2 = []

#calculates bayesian model for all adjacent genes, then plots it and writes it to a file
with open('finalout.txt', 'w') as finalout:	
	pastgene = ''
	for gene in slgenes:
		if pastgene != '':
			finalout.write(pastgene[0] + '\t' +  gene[0] + '\t' + str(evaluateit(density, density2, distance(pastgene[0], gene[0], genetable), fr)) + '\n')
			info.append(evaluateit(density, density2, distance(pastgene[0], gene[0], genetable), fr))
			info2.append(distance(pastgene[0], gene[0], genetable))
		pastgene = gene

axes = plt.gca()
axes.set_xlim([-400,1000])
plt.scatter(info2, info)
plt.xlabel('Gene Distance')
plt.ylabel('Posterior Probability')
plt.savefig('likelyhood.png')

sbpositive = []
sbnegative = []

for op in operonsf:
	prevg = ''
	gio = op[5].split(',')
	for g in gio:
		if prevg != '':
			sbpositive.append(evaluateit(density, density2, distance(g, prevg, genetable), fr))
			
		for pos,gene in enumerate(slgenes):
			if g in gene:
				if slgenes[pos-1][0] not in gio:
					if distance(g, slgenes[pos-1][0], genetable) < 50000:
						sbnegative.append(evaluateit(density, density2, distance(g, slgenes[pos-1][0], genetable), fr))
				if slgenes[pos+1][0] not in gio:
					if distance(g, slgenes[pos+1][0], genetable) < 50000:
						sbnegative.append(evaluateit(density, density2, distance(g, slgenes[pos+1][0], genetable), fr))
		prevg = g

sbpositive.sort()
sbnegative.sort()

spa = []
sea = []
aca = []
scl = []

for thr in range(21):
	tpc = 0
	tnc = 0
	fpc = 0
	fnc = 0
	scl.append(thr/20.0)
	for x in sbpositive:
		if x > thr/20.0:
			tpc += 1
		else:
			fnc += 1
			
	for x in sbnegative:
		if x > thr/20.0:
			fpc += 1
		else:
			tnc += 1
			
	spa.append(tnc/(fpc+tnc))
	sea.append(tpc/(fnc+tpc))
	aca.append((tpc+tnc)/(tpc+tnc+fpc+fnc))
	
plt.close('all')
plt.scatter(scl, spa)
plt.scatter(scl, sea)
plt.scatter(scl, aca)
plt.legend(['Specificity', 'Sensitivity', 'Accuracy'])
plt.xlabel('threshold')

plt.savefig('evaluateresults.png')
		
plt.close('all')

