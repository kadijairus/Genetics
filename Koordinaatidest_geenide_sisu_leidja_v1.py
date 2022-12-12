# Koordinaatide järgi CNV-de sisse jäävate geenide otsija
# Sisend: Tööleht CNVdega, gtf fail geeninimedega, txt fail omim geenidega
# Väljund: Exceli fail valku kodeerivate geenide ja OMIM Morbid tüüp 3 geenide nimede ja arvuga.
# 11.12.2022
# v1
# Kadi Jairus

import os
import pandas as pd
from datetime import datetime
# https://github.com/openvax/gtfparse
from gtfparse import read_gtf
import re

tooleht = input('Sisesta töölehe nimi: ')
if tooleht == '':
    tooleht = 'Tooleht.xlsx'
else:
    tooleht = tooleht
    
print(f'Hakkan otsima geene failis {tooleht}')

# Koordinaatidega fail, millele vaja lisada geeniandmed
#tooleht = 'Tooleht.xlsx'
#https://omim.org/downloads
# Copyright (c) 1966-2022 Johns Hopkins University. Use of this file adheres to the terms specified at https://omim.org/help/agreement
omim = 'mim2gene.txt'
morbid = 'morbidmap.txt'
# https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_42/
# https://www.gencodegenes.org/pages/data_format.html
gtf_fail = 'gencode.v42.chr_patch_hapl_scaff.basic.annotation.gtf'

print("Programm alustas tööd, palun oota...")
df = read_gtf(gtf_fail)

# Index(['seqname', 'source', 'feature', 'start', 'end', 'score', 'strand',
#       'frame', 'gene_id', 'gene_type', 'gene_name', 'level', 'tag',
#       'transcript_id', 'transcript_type', 'transcript_name',
#       'transcript_support_level', 'havana_transcript', 'exon_number',
#       'exon_id', 'hgnc_id', 'havana_gene', 'ont', 'protein_id', 'ccdsid',
#       'artif_dupl'],
#       dtype='object')
#print(df)

# filter DataFrame to gene entries
df_genes = df[df["feature"] == "gene"]

#print('Sellised tulbanimed')
#print(df_genes.keys())
#print('Selline geenide tabel')
#print(df_genes)
#df_genes_chr21 = df_genes[df_genes["seqname"] == "chr21"]
#df_genes_chr21_protein = df_genes_chr21[df_genes_chr21["gene_type"] == "protein_coding"]
#print(df_genes_chr21_protein)
#print(df_genes_chr21_protein.keys())
# Key-value pairs in 9th column: gene_id, gene_type: protein_coding, gene_name

try:
#if True:    
    df_morbid = pd.read_csv(morbid, header = 3, sep='\t', lineterminator='\n')
    df_morbid[['Gene','Gene Symbols']]=df_morbid['Gene Symbols'].str.split(',',n=1,expand=True)  
#    print(df_morbid['Gene'])
    df_morbid[['# Phenotype','Phenotype']]=df_morbid['# Phenotype'].str.split('(',n=1,expand=True)
    df_morbid['Phenotype']=df_morbid['Phenotype'].replace(to_replace='\)',value='',regex=True)
    df_morbid['Phenotype'] = df_morbid['Phenotype'].astype(str)
    #    print(df_morbid['Phenotype'])
    omim_genes = df_morbid['Gene'].tolist()
#    morbid_gene_filter = (df_morbid['Phenotype'] > 3)
    morbid_gene_filter = (df_morbid['Phenotype'] == '3')
    df_morbid_genes = df_morbid[morbid_gene_filter]
#    print(df_morbid_genes['Phenotype'])
    morbid_genes = df_morbid_genes['Gene'].tolist()
#    print(omim_genes)
#    print(morbid_genes)
except:
    print('Probleem OMIMidega')
#i = 1/0

try:
#if True:
    df_tl = pd.read_excel(tooleht)
#    print('Selline tööleht:')
#    print(df_tl)
except:
#else:
    print("Tekkis probleem töölehe faili lugemisel ja filtreerimisel")


def geeniotsija (kromosoom,algus,lopp):
    if str(kromosoom) == '':
        return None
    else:
        df_genes_chr = df_genes[df_genes["seqname"] == 'chr' + str(kromosoom)]
        df_proteincoding = df_genes_chr[df_genes_chr["gene_type"] == "protein_coding"]
        geenisiseste_filter = (((algus <= df_proteincoding['start']) & (lopp >= df_proteincoding['end'])) | ((algus <= df_proteincoding['start']) & (lopp <= df_proteincoding['end']) & (lopp > df_proteincoding['start'])) | ((algus >= df_proteincoding['start']) & (algus < df_proteincoding['end']) & (lopp >= df_proteincoding['end'])))
        df_geenid_sees = df_proteincoding[geenisiseste_filter]
        df_geenid_sees.sort_values(by=['start'])
        genelist = df_geenid_sees['gene_name'].tolist()
#        print('CNV sisse jäävad geenid:')
#        print(genelist)
#        print('OMIM Morbid geenid on:')
#        print(morbid_genes)
        set_genelist = set(genelist)
        set_omimlist = set(morbid_genes)
        omimlist = list(set_genelist.intersection(set_omimlist))
#        omimid_sees_filter = df_geenid_sees['gene_name'].isin(morbid_genes)
#        df_omim_geenid_sees = df_proteincoding[omimid_sees_filter]
#        omimlist = df_geenid_sees['gene_name'].tolist()
        print(omimlist)
        return ((", ".join(genelist)) + ': ' + str(len(set(genelist))) + '/' + (", ".join(omimlist)) + ': ' + str(len(set(omimlist))))


try:
#if True:
    df_tl['Geenid'] = df_tl.apply(lambda x: geeniotsija(x.Kr, x.Algus, x.Lõpp), axis=1)
#    print('Selline tabel enne geenide eri lahtritesse jagamist')
#    print(df_tl['Geenid'])
    df_tl[['Geenid','OMIM Morbid']] = df_tl['Geenid'].str.split('/',n=1,expand=True)
    df_tl[['Geenid','Geenid nr']] = df_tl['Geenid'].str.split(':',n=1,expand=True)
    df_tl[['OMIM Morbid','OMIM nr']] = df_tl['OMIM Morbid'].str.split(':',n=1,expand=True)
except:
    print("Tekkis probleem geeniotsijaga")
    
#print(df_tl['Geenid'])

try:
    praegune_aeg = datetime.now().strftime("%Y-%m-%d_%H-%M")
    tooleht = f"Geeninimedega_{praegune_aeg}_{tooleht}"
    df_tl.to_excel(tooleht,index=False)
    print('Geeninimedega tööleht on tehtud')
except:
    print("Tekkis probleem Exceli faili uuendamisel")