# Carico file
rm(list=ls())
cat("\14")
specie="cow"

#library
library(refGenome)

gtf_path="../Bos_taurus.UMD3.1.88.gtf"
snp_path="snp.txt"

gtf_file="../Bos_taurus.UMD3.1.88.gtf"
snp_file="snp.txt"

#Function
search_genes <- function(gtf_file,snp_file,method=c("gene","exon"),window=250000){
  
  #check method
  method <- match.arg(method)
  message(paste("You are using the method:", method))
  
  #TODO print parametri, file usati
  
  #load file
  #TODO lettura dei path
  if (!(file.exists(gtf_file))){
    stop(paste("file", gtf_file, "doesn't exists"))
  }

  #check SNP file
  if (file.exists(snp_file)){
    markers=read.table(snp_file,header=T)
    colnames(markers)=c("chr","snp","position","value")
  }else{
    stop(paste("file", snp_file, "doesn't exists"))
  }
  
  #read gtf file
  gtf=ensemblGenome()
  read.gtf(gtf,gtf_file)

  #Creating gene data frame
  if (method=="gene"){
    genes=data.frame(gtf@ev$genes)
    gene=genes[,c("seqid","gene_id","gene_name","start","end")]
    colnames(gene)=c("chr","ID","gene_name","start","end")
  }else {
    entire=data.frame(gtf@ev$gtf)
    exon=entire[entire$feature=="exon",]
    gene=exon[,c("seqid","exon_id","gene_name","start","end")]
    colnames(gene)=c("chr","ID","gene_name","start","end")
  }
  

  #Search genes
  gene$lower=gene$start-window
  gene$upper=gene$end+window
  chromosomes=markers$chr
  for (chr in unique(chromosomes)) { # chr in 1:ncrom usando specie
    print(chr)
    temp_gene=gene[gene$chr==chr,]
    temp_markers=markers[markers$chr==chr,]
    temp_gene$ok=F
    posizione=temp_markers[,c("snp","position")]
    if (nrow(temp_markers) > 0) {
      for (k in 1:nrow(posizione)) {
        for (n in 1:nrow(temp_gene)) {
          if (posizione[k,2] >= temp_gene[n,"lower"] & posizione[k,2] <= temp_gene[n,"upper"]) {
            temp_gene[n,"ok"]=T
            temp_gene[n,"snp"]=posizione[k,"snp"]
            temp_gene[n,"snp_position"]=posizione[k,"position"]
          }
        }
      }
    temp_gene=temp_gene[temp_gene$ok,]
    if (chr == chromosomes[1]) {
      final=temp_gene
    }else { 
      final=rbind(final,temp_gene)
      }
    }
  }  
  
  final=final[,c(1:7,9:10)]
  return(final)
}


##### CONTROLLO SE FUNZIONA IL NOME DELLO SNP #####

finale_genes=search_genes(gtf_file = gtf_path ,snp_file = snp_path ,method="gene")
head(finale_genes)

for (i in 1:nrow(finale_genes)) {
  if (finale_genes[i,"snp_position"]<=finale_genes[i,"upper"] & finale_genes[i,"snp_position"]>=finale_genes[i,"lower"]) {
    finale_genes$controllo=T
  }
}




