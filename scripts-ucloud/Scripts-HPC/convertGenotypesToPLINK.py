import sys
from collections import Counter

#infile = "input/139_ouh_june_2017/139_ouh_female_hisp_brca1_onco_geno.txt"
#infile = "input/139_ouh_june_2017/139_ouh_female_hisp_brca2_onco_geno.txt"
# infile = "input/139_ouh_june_2017/139_ouh_brca1_onco_geno.txt"
infile = "input/139_ouh_june_2017/139_ouh_brca2_onco_geno.txt"
outfile = infile.split(".txt")[0]+"_plink_format.txt"

outTable = []

# loadAll = False
loadAll = True


# Count ncols in each row
# with open(infile, "r") as f:
# 	for line in f:
# 		print(len(line.split("\t")))
# 	sys.exit()


def nextCol(start, stop):
	for c in range(start, stop):
		with open(infile, "r") as f:
			col = [x.strip().split("\t")[c] for x in f.readlines()]
			yield c, col

with open(infile, "r") as f:
	first_line = f.readline().strip().split("\t")
	ncol = len(first_line)
	print(ncol)
	
	if loadAll:
		rows = [line.strip().split("\t") for line in f]
		cols = [list(col) for col in zip(*rows)]
	
	print("Finished reading input data")

# cols = [[] for i in range(ncol)]
# with open(infile, "r") as f:
# 	for line in f:
# 		for i, e in enumerate(line.split("\t")):
# 			cols[i].append(e)


# Patient IDs - Remove the "corner" header SNP, as it's already in first_line
header_col = [x.strip().split("\t")[0] for x in open(infile).readlines()][1:]
print(len(header_col))



#for c, col in nextCol(11443,ncol):
#for c in range(17,18):
#for c in range(1,10000):
#for c in range(17,21):
#for c in range(1,4):
for c in range(1, ncol):
	#print(c)
	if c % 500 == 0:
		print("Progress:", c, "of", ncol, "columns")

	
	#col = [x.split("\t")[c] for x in f.readlines()]
	#col = cols[c]
	
	#idx = col[0]
	#genotypes = col[1:]
	#print(idx)
	#print(genotypes)
	
	genotypes = cols[c]
	
	baseCounts = Counter()
	for geno in genotypes:
		for base in geno:
			if base == "-":
				continue
			baseCounts[base] += 1
	
	#print(baseCounts)
	
	plink_col = []
	if len(baseCounts) == 1:
		for i in range(len(genotypes)):
			plink_col.append(0)


	if len(baseCounts) == 2:
		A_allele = ("", 0)
		B_allele = ("", 0)
		
		for base in baseCounts:
			if baseCounts[base] > A_allele[1] or A_allele[1] == 0:
				A_allele = (base, baseCounts[base])
			if baseCounts[base] < B_allele[1] or B_allele[1] == 0:
				B_allele = (base, baseCounts[base])
		
		#print(A_allele)
		#print(B_allele)
		
		homoA = A_allele[0]+A_allele[0]
		hetAB = A_allele[0]+B_allele[0]
		hetBA = B_allele[0]+A_allele[0]
		homoB = B_allele[0]+B_allele[0]
		
		for geno in genotypes:
			if geno == homoA:
				plink_col.append(0)
			elif geno == hetAB:
				plink_col.append(1)
			elif geno == hetBA:
				plink_col.append(1)
			elif geno == homoB:
				plink_col.append(2)
			else:
				plink_col.append(-1)
				#print(c)
				#print(geno)
		

	if len(baseCounts) > 2:
		print("CRAP!")
		print(c)
		print(col)
		print(baseCounts)
		sys.exit()
	
	# print(plink_col)
	
	outTable.append(plink_col)

with open(outfile, "w") as f:
	f.write('\t'.join(str(e) for e in first_line)+'\n')
	for row in zip(header_col, *outTable):
	#for row in zip(cols[0], *outTable):
		# print row
		f.write('\t'.join(str(r) for r in row)+'\n')