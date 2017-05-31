from Bio import SeqIO

repe = []
repa = []

x = -1

for record in SeqIO.parse('ecoli.gbff', 'gb'):

	#each dataset has multiple features
	for feature in record.features:
	
		#gets the taxid, and replicon type
		if feature.type == 'source':
			x += 1
		
		#get information about each gene
		if feature.type == 'CDS':
			j = feature.qualifiers.get('protein_id', '')
			if j != '':
				repe.append([j])

x = -1

for record in SeqIO.parse('atume.gbff', 'gb'):

	#each dataset has multiple features
	for feature in record.features:
	
		#gets the taxid, and replicon type
		if feature.type == 'source':
			x += 1
			repa.append([])
		
		#get information about each gene
		if feature.type == 'CDS':
			j = feature.qualifiers.get('protein_id', '')
			if j != '':
				repa[x].append([j])

pairdict = {}

with open('finalpairs.txt', 'r') as fp:
	for line in fp:
		l = line.split('\t')
		pairdict[l[1].lstrip().rstrip()] = l[0].lstrip().rstrip()
		
final = []
		
for atgene in pairdict.keys():
	for replicon in repa:
		for c,x in enumerate(replicon):
			print('hi')
			if atgene in x:
				for i in range(c-5, c+5):
					print('hi')
					if i >= 0 and i < len(replicon) and i != c:
						print('2')
						if replicon[i] in pairdict.keys():
							for q in range(len(repe)):
								if pairdict(replicon[i]) in repe[q]:
									for z in range(q-5, c+5):
										if z >= 0 and z < len(replicon) and z != c:
											if pairdict[atgene] in repe[z]:
												final.append([atgene,pairdict[atgene]])
												
print(final)
print(len(final))
		
