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
###
library(tidyverse)
library(coin)
###Chi-square test for metabolic pathway pairwise between tumor and normal tissues
for (l in 1:13){
  dat<-vector(mode = "list")
  p_value<-vector()
  or <- vector()
  tcga<-tcga_list[l]
  #cancer<-read.table(file = paste("H:/tcga/model/",tcga,".txt",sep = ""),sep = "\t",header = T,check.names = F)
  path_pairs<-read.table(file = paste("D:/tcga/model/reconstruct_tpm/tpm/tpm_dif/",tcga,"paths_dif.txt",sep = ""),sep = "\t",header = T,check.names = F)
  typ = substring(colnames(path_pairs),14,16)
  #print(table(typ))
  wang <- 0.1
  for (j in 1:3486) {
    x = 0
    dif <- as.numeric(path_pairs[j,])
    type <- typ
    if(sum(dif <= -wang | dif >= wang) >= length(dif)*0.8){
      dif <- dif[dif <= -wang | dif >= wang]
      type <- typ[dif <= -wang | dif >= wang]
      print("delete");print(j)
    }
    
    dif = ifelse(dif > 0, 0, 1)
      
      a<-dif[type == "01A"] %>% .[. == 0] %>% length()
      b<-dif[type == "11A"] %>% .[. == 0] %>% length()
      c<-dif[type == "01A"] %>% .[. == 1] %>% length()
      d<-dif[type == "11A"] %>% .[. == 1] %>% length()
      
      ###Fisher exact test for specific metabolic pathway pairwise comparison
      q<-min(((a+c)*(a+b)/(a+b+c+d)),((b+d)*(a+b)/(a+b+c+d)),((a+b)*(d+b)/(a+b+c+d)),((d+c)*(d+b)/(a+b+c+d)))
      if(q < 5){print(c(i,j))
        x<-x+1
        dat[[x]]<-c(i,j)
        chisq1<-fisher.test(matrix(c(a,b,c,d),nrow = 2))$p.value
        #print(matrix(c(a,b,c,d),nrow = 2))
        chisq2<-fisher.test(matrix(c(c,d,a,b),nrow = 2))$p.value 
      }
      else{
        chisq1<-chisq_test(as.table(matrix(c(a,b,c,d),nrow = 2)),distribution = "exact") %>% pvalue() %>% as.numeric()
        #print(matrix(c(a,b,c,d),nrow = 2))
        chisq2<-chisq_test(as.table(matrix(c(c,d,a,b),nrow = 2)),distribution = "exact") %>% pvalue() %>% as.numeric()
      }
      or1 = (a*d)/(b*c)
      or2 = (b*c)/(a*d)
      p_value<-append(p_value,c(chisq1,chisq2))
      or <- append(or,c(or1,or2))
  }   
  #lst<-do.call("rbind",lst)
  chiq_adjust1<-p.adjust(p_value[seq(1,6971,2)],method = "fdr")
  chiq_adjust2<-p.adjust(p_value[seq(2,6972,2)],method = "fdr")
  write.table(cbind(p_value[seq(1,6971,2)],p_value[seq(2,6972,2)],chiq_adjust1,chiq_adjust2,or[seq(1,6971,2)],or[seq(2,6972,2)]),file = paste("D:/tcga/model/reconstruct_tpm/tpm/tpm_test/",tcga,"path_pairs.txt",sep = ""),sep = "\t",col.names = F,row.names = F)
  dat<-do.call("rbind",dat)
  if(length(dat > 0)){write.table(dat,file = paste("D:/tcga/model/reconstruct_tpm/tpm/",tcga,"path_fisher.txt",sep = ""),sep = "\t",col.names = F,row.names = F)
  }}
####Cramer's coefficient
#tcga_list<- c( "TCGA-BRCA", "TCGA-COAD", "TCGA-HNSC", "TCGA-KICH", "TCGA-KIRC", "TCGA-KIRP", "TCGA-LIHC","TCGA-LUAD", "TCGA-LUSC", "TCGA-PRAD", "TCGA-STAD", "TCGA-THCA", "TCGA-UCEC")
x = 0
datv<-vector(mode = "list")
for (l in 1:13){

  tcga<-tcga_list[l]
  #cancer<-read.table(file = paste("H:/tcga/model/",tcga,".txt",sep = ""),sep = "\t",header = T,check.names = F)
  path_pairs<-read.table(file = paste("D:/tcga/model/reconstruct_tpm/tpm/tpm_dif/",tcga,"paths_dif.txt",sep = ""),sep = "\t",header = T,check.names = F)
  
  typ = substring(colnames(path_pairs),14,16)
  #print(table(typ))
  wang <- 0.1
  
  v<-vector()
  for (j in 1:3486) {
    
    
    dif <- as.numeric(path_pairs[j,])
    type <- typ
    if(sum(dif <= -wang | dif >= wang) >= length(dif)*0.8){
      dif <- dif[dif <= -wang | dif >= wang]
      type <- typ[dif <= -wang | dif >= wang]
      print("delete");x<-x+1;print(x)
    }
    
    dif = ifelse(dif > 0, 0, 1)
    
    a<-dif[type == "01A"] %>% .[. == 0] %>% length()
    b<-dif[type == "11A"] %>% .[. == 0] %>% length()
    c<-dif[type == "01A"] %>% .[. == 1] %>% length()
    d<-dif[type == "11A"] %>% .[. == 1] %>% length()
    v <-append(v, assocstats(matrix(c(a,b,c,d),nrow = 2))$cramer)
  }
  datv[[l]]<-v
  }
