import pandas as pd
import numpy as np
import os

path = 'D:/OneDrive - Systems Oncology, LLC/Bioinformatics/'
pathname = os.fspath(path)

filename = 'pivoted_CERES_FC_w_DepMap_ID.csv'

filepathname = pathname + filename

df_wideT = pd.read_csv(filepathname, index_col=0) #wide format, cell lines index, gene symbols as columns
df_wide = df_wideT.T 
df1 = pd.read_csv('genome_annotations_biomart_updated.csv')
dgd = pd.read_csv('dgd_Hsa_all_v71.tsv', sep='\t', index_col='Name')

print(df3.columns.values)
#load Biomart file with updated gene names
df2 = pd.read_csv('mart_export_updated.txt', delimiter='\t')
df2.head()

#replace with new gene names
df3 = pd.merge(df1, df2, left_on= "Stable_ID", right_on="StableID", how = 'left')
df3['New_Gene'].fillna(df3['Gene'], inplace= True)
df3 = df3.drop(columns=['StableID','Gene'])
df3.head()

df3.rename(columns={'New_Gene' : 'Gene'}, inplace=True)
df3.reindex(columns= ['Gene','Stable_ID','Band','GCpercent','Start','End','Strand','Chromosome'])

#df3.loc[df3.New_Gene == 'A1BG'] #find index of particular gene

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
            chromosome = df3[df3['Gene'] == gene].Chromosome.values[0]
            df5 = df3[df3['Chromosome'] == chromosome].sort_values(by='End').reset_index()
            pos = df5[df5['Gene'] == gene].index[0]
            locus_genes = df5.loc[pos-rank_threshold:pos+rank_threshold, :].Gene.tolist()
            locus_genes.remove(gene)

        a[gene] = df_wide[gene] - (df_wide[[i for i in locus_genes if i in df_wide]].median(axis=1)/2)

with open('missing_genes.txt', 'w') as filehandle:
    for listitem in missing_genes:
        filehandle.write('%s\n' % listitem)

b = a.corr()
b.to_hdf('corrected_coessentiality_matrix.h5', key='stage', mode='w')
print('Done')