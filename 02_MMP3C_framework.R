###13 cancer types of gene expression filename
tcga_list<- c( "TCGA-BRCA", "TCGA-COAD", "TCGA-HNSC", "TCGA-KICH", "TCGA-KIRC", "TCGA-KIRP", "TCGA-LIHC","TCGA-LUAD", "TCGA-LUSC", "TCGA-PRAD", "TCGA-STAD", "TCGA-THCA", "TCGA-UCEC")
files<-list.files(path = "D:/tcga/fpkm/",pattern = ".gz$")
files<-files[match(tcga_list,substring(files,1,9))]
###gene annotation files
anno<-read.table(file = "D:/tcga/anno_file/gencode.v22.annotation.gene.probeMap",sep = "\t",header = T)
###gene symbol annotation
dat<-vector(mode = "list")
for (i in 1:13) {
  cancer<-read.table(file = paste("D:/tcga/fpkm/",files[i],sep = ""),sep = "\t",header = T,check.names = F)
  cancer<-name_ver(cancer)
  print(table(substring(colnames(cancer),14,16)))
  cancer<-cancer[,substring(colnames(cancer),14,16) == "01A"]
  cancer <- 2^cancer-1
  cancer <- apply(cancer,2,fpkmToTpm)
  cancer <- log2(cancer + 1)
  print(head(cancer[1:3,1:4]))
  dat[[i]]<-cancer
}
###load KEGG pathway annotation data and gene weight data computed by Page rank method
load(file = "H:/tcga/path/keggpath.Rdata")
load(file = "D:/tcga/bc/me_weight.Rdata")

hsa_gene<-vector()
for (i in 1:84){
  hsa_name <- get(keggpath[[i]][[1]]$ENTRY)
  hsa_gene <- append(hsa_gene,hsa_name$V1)
}

###
pas<-vector(mode = "list")
for (i in 1:13) {
  hsa_mat <- vector(mode = "list")
  cancer<-dat[[i]]
  print(dim(cancer))
  print(table(hsa_gene %in% row.names(cancer)))
  print(max(cancer[match(hsa_gene,row.names(cancer)),]))
  print(min(cancer[match(hsa_gene,row.names(cancer)),]))
  for (j in 1:84) {
    hsa_name <- get(keggpath[[j]][[1]]$ENTRY)
    can <- cancer[match(hsa_name$V1,row.names(cancer)),]
    score_mat <- t(can) %*% hsa_name$V2
    hsa_mat[[j]]<-score_mat[,1]
  }
  hsa_mat <- do.call("rbind",hsa_mat)
  #write.table(hsa_mat,file = paste("D:/tcga/model/reconstruct_tpm/tpm/tpm_pas/",tcga_list[i],"pas.txt",sep = ""),col.names = T,row.names = F,sep = "\t",quote = F)
  pas[[i]]<-hsa_mat
  }
### Normlization of metabolic pathway activity score
pas_norm <- vector(mode = "list")
for (i in 1:13) {
  cancer<-pas[[i]]
  print(max(cancer));print(min(cancer))
  can <- apply(cancer, 1, FUN = function(x){
    lu1 <- max(x)
    lu2 <- min(x)
    lu3 <- (x - lu2)/(lu1 - lu2)
    return(lu3)
  })
  can <- t(can)
  pas_norm[[i]]<-can
  #write.table(can,file = paste("D:/tcga/model/reconstruct_tpm/tpm/norm/",tcga_list[i],"pas_norm.txt",sep = ""),col.names = T,row.names = F,sep = "\t",quote = F)
}
###the comparison between two metabolic pathways
dif_lst<-vector(mode = "list")
for (l in 1:13){
  wang<-vector(mode = "list")
  tcga<-tcga_list[l]
  path_pairs<-pas_norm[[l]]
  type = substring(colnames(path_pairs),14,16)
  print(table(type))
  x = 0
  for (i in 1:83) {
    for (j in (i+1):84) {
      dif = path_pairs[i,]-path_pairs[j,]
      x = x + 1
      wang[[x]]<-dif
    }}
  wang<-do.call("rbind",wang)
  colnames(wang) <- colnames(path_pairs)
  dif_lst[[l]]<-wang
  #write.table(dat,file = paste("D:/tcga/model/reconstruct_tpm/tpm/tpm_dif/",tcga,"paths_dif.txt",sep = ""),sep = "\t",col.names = T,row.names = F,quote = F)
}
