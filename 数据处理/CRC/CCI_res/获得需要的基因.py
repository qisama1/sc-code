gene1 = set()

for i in c1.index:
    idx = i
    if (len(idx.split('_')) > 1):
        idx = idx.split('_')[0]
        gene_c = i.split('_')[1]
        gene1.add(gene_c)
    gene1.add(idx.split('|')[0])
    gene1.add(idx.split('|')[1])

gene2 = set()

for i in c2.index:
    idx = i
    if (len(idx.split('_')) > 1):
        idx = idx.split('_')[0]
        gene_c = idx.split('_')[1]
        gene2.add(gene_c)
    gene2.add(idx.split('|')[0])
    gene2.add(idx.split('|')[1])
