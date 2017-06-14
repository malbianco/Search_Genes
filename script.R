# Carico file
rm(list=ls())
cat("\14")
specie="cow"

#library
library(refGenome)


gtf="Bos_taurus.UMD3.1.88.gtf"
snp_file="snp.txt"
method="gene"


#Function
search_genes <- function(gtf_file,snp_file,method=c("gene","exon"),window=250000){
  
  #check method
  method <- match.arg(method)
  message(paste("You are using the method:", method))
  
  #load file
  #TODO lettura dei path
  gtf=read.table(paste0("../",gtf),fill=T)
  markers=read.table(snp_file,header=T)
  colnames(markers)=c("chr","snp","position","value")

  #Per ora lo lasciamo
  # if (specie=="cow") ncrom="29"
  # if (specie=="sheep") ncrom="26"
  
  ## start function
  if (method == "gene") {
    #print(paste("searching:",method))
    gene=subset(gtf,V3=="gene",select=c(1,10,3,4,5))
    colnames(gene)=c("chr","name","type","start","end")
    gene$lower=gene$start-window
    gene$upper=gene$end+window
    for (chr in markers$chr) { # chr in 1:ncrom usando specie
      temp_gene=gene[gene$chr==chr,]
      temp_markers=markers[markers$chr==chr,]
      temp_gene$ok=F
      if (nrow(temp_markers) > 0) {
        for (k in 1:nrow(temp_markers)) {
          posizione=temp_markers[k,"position"]
          for (n in 1:nrow(temp_gene)) {
            if (posizione >= temp_gene[n,"lower"] & posizione <= temp_gene[n,"upper"])
              temp_gene[n,"ok"]=T
        }
      }
      temp_gene=temp_gene[temp_gene$ok,]
      if (chr == 1) genes=temp_gene
      else genes=rbind(genes,temp_gene)
   }
  }
 } 
  else {
    #print(paste("searching:",method))
    gene=subset(gtf,V3=="exon",select=c(1,10,3,4,5))
    colnames(gene)=c("chr","name","type","start","end")
    gene$lower=gene$start-window
    gene$upper=gene$end+window
    for (chr in markers$chr) { # chr in 1:ncrom usando specie
      temp_gene=gene[gene$chr==chr,]
      temp_markers=markers[markers$chr==chr,]
      temp_gene$ok=F
      if (nrow(temp_markers) > 0) {
        for (k in 1:nrow(temp_markers)) {
          posizione=temp_markers[k,"position"]
          for (n in 1:nrow(temp_gene)) {
            if (posizione >= temp_gene[n,"lower"] & posizione <= temp_gene[n,"upper"])
              temp_gene[n,"ok"]=T
          }
        }
        temp_gene=temp_gene[temp_gene$ok,]
        if (chr == 1) genes=temp_gene
        else genes=rbind(genes,temp_gene)
      }
    }
  } 
genes=genes[,1:5]
return(genes)
}



genes=search_genes(gtf,snp_file,method="gene")
