import pandas as pd

df = pd.read_csv("./desq2_Sanchez.csv", sep="\t")
# groups = df.groupby(['gene_name', 'baseMean', 'log2FoldChange', 'lfcSE', 'stat', 'pvalue', 'padj']) 

df2 = pd.read_csv("./tsses.gtf", sep="\t" ,header= None)

df2_names =  []
gene_used = []
gene_total = []
gene_repetidos = []

cont = 0

print( df)

for i in df2[8]:
    df2_names.append(i.split(";")[3].split(" ")[2].split('"')[1])

for x in df['gene_name']:
    if x in gene_total:
        gene_repetidos.append(x)
        
    else:
        gene_total.append(x)
    if x in gene_used:
        continue
    else:
        gene_used.append(x)
        if x in df2_names:
            cont = cont +1

print("Contador total mixto entre ambos: ",cont, "CAntidad de elementos en deseq2: ", len(df['gene_name']),"Genes iguales en deseq2: ", gene_repetidos)

