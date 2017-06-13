# Carico file
rm(list=ls())
cat("\14")
specie="cow"
window=250000
gtf="Bos_taurus.UMD3.1.88.gtf"
snp="snp.txt"
method="gene"

genes=fai(gtf,snp,method)

## ho solo un dubbio sulla finestra perchè nei primi geni della lista
## crea un valore negativo per lower

fai=function(gtf,snp,type=c("gene","exon"))
{
  gt=read.table(paste0("../",gtf),fill=T)
  markers=read.table(snp,header=T)
  colnames(markers)=c("chr","snp","position","value")
  if (specie=="cow") ncrom="29"
  if (specie=="sheep") ncrom="26"
  ## start function
  if (method == "gene") {
    print(paste("searching:",method))
    gene=subset(gt,V3=="gene",select=c(1,10,3,4,5))
    colnames(gene)=c("chr","name","type","start","end")
    gene$lower=gene$start-window
    gene$upper=gene$end+window
    for (chr in 1:ncrom) { # chr in 1:ncrom usando specie
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
    print(paste("searching:",method))
    gene=subset(gt,V3=="gene",select=c(1,10,3,4,5))
    colnames(gene)=c("chr","name","type","start","end")
    gene$lower=gene$start-window
    gene$upper=gene$end+window
    for (chr in 1:ncrom) { # chr in 1:ncrom usando specie
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


