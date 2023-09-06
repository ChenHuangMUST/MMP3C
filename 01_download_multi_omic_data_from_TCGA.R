library(downloader)
tcga_list<- c( "TCGA-BRCA", "TCGA-COAD", "TCGA-HNSC", "TCGA-KICH", "TCGA-KIRC", "TCGA-KIRP", "TCGA-LIHC","TCGA-LUAD", "TCGA-LUSC", "TCGA-PRAD", "TCGA-STAD", "TCGA-THCA", "TCGA-UCEC")
#https://gdc-hub.s3.us-east-1.amazonaws.com/download/TCGA-BRCA.survival.tsv
for (i in 1:13){
  tcga<-tcga_list[i]
  download(paste("https://gdc-hub.s3.us-east-1.amazonaws.com/download/",tcga,".methylation450.tsv.gz",sep = ""),
           paste(tcga,".masked_cnv.tsv.gz",sep = ""), mode = "wb")
}
#https://gdc-hub.s3.us-east-1.amazonaws.com/download/TCGA-LIHC.htseq_fpkm.tsv.gz
###cnv download
#https://gdc-hub.s3.us-east-1.amazonaws.com/download/TCGA-PRAD.masked_cnv.tsv.gz
#https://gdc-hub.s3.us-east-1.amazonaws.com/download/TCGA-LUAD.mutect2_snv.tsv.gz
#https://gdc-hub.s3.us-east-1.amazonaws.com/download/TCGA-LIHC.methylation450.tsv.gz
