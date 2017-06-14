# Carico file
rm(list=ls())
cat("\14")
specie="cow"

#library
library(refGenome)


gtf_file="Bos_taurus.UMD3.1.88.gtf"
snp_file="snp.txt"
window=250000

#Function
search_genes <- function(gtf_file,snp_file,method=c("gene","exon"),window=250000){
  
  #check method
  method <- match.arg(method)
  message(paste("You are using the method:", method))
  
  #load file
  #TODO lettura dei path
  gtf=ensemblGenome()
  read.gtf(gtf,paste0("../",gtf_file))
  genes=data.frame(gtf@ev$genes)
  entire=data.frame(gtf@ev$gtf)
  exon=entire[entire$feature=="exon",]
  markers=read.table(snp_file,header=T)
  colnames(markers)=c("chr","snp","position","value")

  #Per ora lo lasciamo
  # if (specie=="cow") ncrom="29"
  # if (specie=="sheep") ncrom="26"
  
  ## start function
  if (method == "gene") {
    #print(paste("searching:",method))
    gene=genes[,c("seqid","gene_id","gene_name","start","end")]
    colnames(gene)=c("chr","gene_id","gene_name","start","end")
    gene$lower=gene$start-window
    gene$upper=gene$end+window
    chromosomes=markers$chr
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
      if (chr == chromosomes[1]) final=temp_gene
      else final=rbind(final,temp_gene)
   }
  }
 } 
  else {
    #print(paste("searching:",method))
    gene=exon[,c("seqid","exon_id","gene_name","start","end")]
    colnames(gene)=c("chr","exon_id","gene_name","start","end")
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
        if (chr == chromosomes[1]) final=temp_gene
        else final=rbind(final,temp_gene)
      }
    }
  } 
final=final[,1:5]
return(final)
}



genes=search_genes(gtf_file,snp_file,method="gene")




