# Find SNP position by mergin on variants.tsv from Neale Lab

import pandas as pd
import os

def main():


    variants = pd.read_csv('/Users/jballard/Documents/FMR_multi/gene_enrichment/data/variants.tsv',sep='\t',low_memory=False)
    variants = variants[["chr","pos","rsid"]]

    fprefix = 'coronary_3fr_082621'


    for i in range(1,10):
        fname = '/Users/jballard/Documents/FMR_multi/gene_enrichment/data/'+fprefix+'_cpt'+str(i)+'.csv'

        # First check if file is empty
        if os.stat(fname).st_size == 0:
            continue

        cpt_snplist = pd.read_csv(fname,header=None)
        cpt_snplist.columns=['rsid']
        cpt_snpjoin = pd.merge(variants,cpt_snplist)
        cpt_snpjoin.to_csv('/Users/jballard/Documents/FMR_multi/gene_enrichment/data/'+fprefix+'_cpt'+str(i)+'_snppos.csv',sep=',',index=False)


if __name__=="__main__":
    main()
