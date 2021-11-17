# Importing list of rsids

library(ACME)

for (c in 1:9) {
  
  tname = 'coronary_3fr_082621'
  
  fname <- paste('/Users/jballard/Documents/FMR_multi/gene_enrichment/data/',tname,'_cpt',c,'_snppos.csv',sep='')

  if (file.exists(fname)) {
  
    rs <- read.csv(fname)
  
    genes <- c()
    for (i in 1:dim(rs)[1]) {
      chrom <- rs$chr[i]
      chromstr <- paste('chr',chrom,sep="")
      pos <- rs$pos[i]
      
      fcg <- findClosestGene(chromstr,pos,genome = "hg19")
      gname <- fcg$geneName[1] # Just taking the first gene name that comes up - some are repeats anyway
      genes <- c(genes, gname)
    }
    
    # Only keeping unique gene names
    genes <- unique(genes)
    
    genetb <- as.data.frame(genes)
    write.table(genetb, paste('/Users/jballard/Documents/FMR_multi/gene_enrichment/data/',tname,'_cpt',c,'_genes.csv',sep=''),quote=FALSE,row.names=FALSE)
  }
}
