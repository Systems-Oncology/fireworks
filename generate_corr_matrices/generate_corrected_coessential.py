import pandas as pd
import numpy as np
import os

path = '~/Systems Oncology, LLC/Mine - Documents/Projects/EXP_002_Differential_Fitness_Oncogenic_Driver_Mut/Design/'
pathname = os.fspath(path)

filename = 'pivoted_CERES_FC_w_DepMap_ID.csv'

filepathname = pathname + filename

df_wideT = pd.read_csv(filepathname, index_col=0) #wide format, cell lines index, gene symbols as columns
df_wide = df_wideT.T 
df1 = pd.read_csv('genome_annotations_biomart_updated.csv')
dgd = pd.read_csv('dgd_Hsa_all_v71.tsv', sep='\t', index_col='Name')

rank_threshold = 20 #on each side
dup_genes = dgd.index.tolist()

missing_genes = []

a = pd.DataFrame()
for gene in df_wide.columns.tolist():
    if gene in dup_genes:
        a[gene] = df_wide[gene] #just copy the original essentiality score if a duplicate gene
    else:
        try: 
            chromosome = df1[df1['Gene'] == gene].Chromosome.values[0]
        except: 
            print(gene)
            missing_genes.append(gene)
        else:
            df2 = df1[df1['Chromosome'] == chromosome].sort_values(by='End').reset_index()
            pos = df2[df2['Gene'] == gene].index[0]
            locus_genes = df2.loc[pos-rank_threshold:pos+rank_threshold, :].Gene.tolist()
            locus_genes.remove(gene)

            a[gene] = df_wide[gene] - (df_wide[[i for i in locus_genes if i in df_wide]].median(axis=1)/2)

 with open('missing_genes.txt', 'w') as filehandle:
    for listitem in missing_genes:
        filehandle.write('%s\n' % listitem)

b = a.corr()
b.to_hdf('corrected_coessentiality_matrix.h5', key='stage', mode='w')
print('Done')