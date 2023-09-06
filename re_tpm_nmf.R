dir.create("D:/tcga/model/reconstruct_tpm/tpm/nmf")
load(file = "D:/tcga/model/correlation2.Rdata")
rm(list = setdiff(ls(),c("name_ver","fpkmToTpm")))
tcga_list<- c( "TCGA-BRCA", "TCGA-COAD", "TCGA-HNSC", "TCGA-KICH", "TCGA-KIRC", "TCGA-KIRP", "TCGA-LIHC","TCGA-LUAD", "TCGA-LUSC", "TCGA-PRAD", "TCGA-STAD", "TCGA-THCA", "TCGA-UCEC")
###
load(file = "D:/tcga/model/correlation2.Rdata")
rm(list = setdiff(ls(),c("name_ver","fpkmToTpm")))
#pancreas<-read.table(file = "H:/tcga/count/TCGA-PAAD.htseq_counts.tsv.gz",sep = "\t",header = T,check.names = F)
files<-list.files(path = "D:/tcga/fpkm/",pattern = ".gz$")
files<-files[match(tcga_list,substring(files,1,9))]

anno<-read.table(file = "D:/tcga/anno_file/gencode.v22.annotation.gene.probeMap",sep = "\t",header = T)
dat<-vector(mode = "list")
for (i in 1:13) {
  cancer<-read.table(file = paste("D:/tcga/fpkm/",files[i],sep = ""),sep = "\t",header = T,check.names = F)
  cancer<-name_ver(cancer)
  print(table(substring(colnames(cancer),14,16)))
  cancer<-cancer[,substring(colnames(cancer),14,16) == "01A"]
  cancer <- 2^cancer-1
  #cancer<-cancer[,colnames(can)]
  cancer <- apply(cancer,2,fpkmToTpm)
  cancer <- log2(cancer + 1)
  print(head(cancer[1:3,1:4]))
  dat[[i]]<-cancer
  #write.table(cancer,file = paste("D:/tcga/model/reconstruct_tpm/tpm_data/",tcga_list[i],"tpm.txt",sep = ""),col.names = T,row.names = T,sep = "\t",quote = F)
}
######
load(file = "H:/tcga/path/keggpath.Rdata")
load(file = "D:/tcga/bc/me_weight.Rdata")

hsa_gene<-vector()

for (i in 1:84){
  #hsa_name<-read.table(
  #  file = paste("H:/tcga/84_path/",keggpath[[i]][[1]]$ENTRY,"/weight.txt",sep = ""),
  #  sep = "\t",header = F)
  #assgin(keggpath[[i]][[1]]$ENTRY,hsa_name)
  hsa_name <- get(keggpath[[i]][[1]]$ENTRY)
  #print(table(hsa_name$V1 %in% row.names(cancer)))
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
    #hsa_mat <- vector(mode = "list")
    hsa_name <- get(keggpath[[j]][[1]]$ENTRY)
    can <- cancer[match(hsa_name$V1,row.names(cancer)),]
    score_mat <- t(can) %*% hsa_name$V2
    hsa_mat[[j]]<-score_mat[,1]
  }
  #hsa_mat[[i]]<-score_mat[,1]
  hsa_mat <- do.call("rbind",hsa_mat)
  #write.table(hsa_mat,file = paste("D:/tcga/model/reconstruct_tpm/tpm/tpm_pas/",tcga_list[i],"pas.txt",sep = ""),col.names = T,row.names = F,sep = "\t",quote = F)
  pas[[i]]<-hsa_mat
  }
### Normlization
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
###dif
####dif
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
######nmf 
load(file = "tpm/reconstruct_tpm_tu.Rdata")

marker <- lapply(can, function(x){
  lu <- log2(x$or1)
  lu[is.infinite(lu)] <- max(lu[!is.infinite(lu)])
  lu <- (lu - mean(lu))/sd(lu)
  index <- (lu > 2 | lu < -2) &  x$v > 0.6
  index <- index & (x$P < 10^(-6))
  return(x$SNP[index])
}) %>% unlist() %>% unique()
marker <- as.numeric(substring(marker,5))
### NMF
library(NMF)

system.time(res <- nmf(nmf_imput[marker,]+1, 2:7, nrun= 20))
res.multi.method1 <- nmf(nmf_imput[marker,]+1, 4,list('brunet', 'lee', 'ns'), seed=123, .options='t')
plot(res.multi.method1)
system.time(res4 <- nmf(nmf_imput[marker,]+1, 4, method = "lee", .options='t'))
system.time(res5 <- nmf(nmf_imput[marker,]+1, 5, method = "lee", .options='t'))
res.multi.method2 <- nmf(nmf_imput[marker,]+1, 5,list('brunet', 'lee', 'ns'), seed=123, .options='t')
plot(res.multi.method2)
####
save.image(file = "tpm/nmf/nmf_1.Rdata")
save.image(res5,file = "tpm/nmf/nmf_2.Rdata")
####OS time survival data
library(NMF)
library(survival)
library(survminer)
sur_data <- vector(mode = "list")
for (i in 1:13) {
  sur<-read.table(file = paste("H:/tcga/survival/",tcga_list[i],".survival.tsv",sep = ""),sep = "\t",header = T,check.names = F)
  sur_data[[i]]<-sur[sur$sample %in% colnames(nmf_imput),]
}
sur_data <- do.call("rbind",sur_data)
####
ac <- data.frame(sample = colnames(nmf_imput), group = predict(res4))
sur_data$group1 <- as.factor(ac$group[match(sur_data$sample,ac$sample)])
sur_data$OS.time <- sur_data$OS.time/365
fit <- survfit(Surv(OS.time,OS) ~ group2,  
               data = sur_data)
fit
ggsurvplot(fit, 
           data = sur_data, 
           conf.int = F,
           pval = TRUE, 
           surv.median.line = "hv", 
           palette = "hue")
#####validation single cancer type
####add single cell type group info
nmf_sur_df <- vector(mode = "list")
sur_data <- vector(mode = "list")
for (i in 1:13) {
  sur<-read.table(file = paste("H:/tcga/survival/",tcga_list[i],".survival.tsv",sep = ""),sep = "\t",header = T,check.names = F)
  sur$group <- i
  sur_data[[i]]<-sur[sur$sample %in% colnames(nmf_imput),]
}
sur_data <- do.call("rbind",sur_data)
##single cancer type sur analysis
lr_p <- vector()
for (i in 1:13) {
  pdf(file = paste("D:/tcga/model/reconstruct_tpm/tpm/nmf/sur5",i,".pdf",sep = ""),width = 12,height = 10)
  dats<-sur_data[sur_data$group == i,]
  
  fit <- survfit(Surv(OS.time,OS) ~ group2,  
                 data = dats)
  #fit
  p<-ggsurvplot(fit, 
                data = dats, 
                conf.int = F,
                pval = TRUE, 
                surv.median.line = "hv", 
                palette = "hue")
  print(p)
  dev.off()
  lr<-survdiff(Surv(OS.time, OS) ~ group2, data = dats)
  lr_p<-append(lr_p,lr$pvalue)
}
###
nmf_sur_df[[2]]<-lr_p
###
nmf_sur_df <- vector(mode = "list")
lr_p <- vector()
  for (i in 6:13) {
    #lr_p <- vector()
    #pdf(file = paste("D:/tcga/model/reconstruct_tpm/tpm/nmf/sur5",i,".pdf",sep = ""),width = 12,height = 10)
    dats<-sur_data[sur_data$group == i,]
    print(table(dats$group2))
    query <- (1:5)[as.numeric(table(dats$group2)) >= 30]
    dats<-dats[dats$group2 %in% query,]
    lr<-survdiff(Surv(OS.time, OS) ~ group2, data = dats)
    lr_p<-append(lr_p,lr$pvalue)
  }
  nmf_sur_df[[2]]<-lr_p
#####
 lr_p<-append(lr_p,NA,after = 3)
#####
