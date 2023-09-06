###The identification of the significant metabolic pathway pairs using MMP3C framework
marker <- lapply(can, function(x){
  lu <- log2(x$or1)
  lu[is.infinite(lu)] <- max(lu[!is.infinite(lu)])
  lu <- (lu - mean(lu))/sd(lu)
  index <- (lu > 2 | lu < -2) &  x$v > 0.6
  index <- index & (x$P < 10^(-6))
  return(x$SNP[index])
}) %>% unlist() %>% unique()
marker <- as.numeric(substring(marker,5))
###NMF clustering
library(NMF)
system.time(res <- nmf(nmf_imput[marker,]+1, 2:7, nrun= 20))
res.multi.method1 <- nmf(nmf_imput[marker,]+1, 4,list('brunet', 'lee', 'ns'), seed=123, .options='t')
plot(res.multi.method1)
system.time(res4 <- nmf(nmf_imput[marker,]+1, 4, method = "lee", .options='t'))
system.time(res5 <- nmf(nmf_imput[marker,]+1, 5, method = "lee", .options='t'))
res.multi.method2 <- nmf(nmf_imput[marker,]+1, 5,list('brunet', 'lee', 'ns'), seed=123, .options='t')
plot(res.multi.method2)
###
save.image(file = "tpm/nmf/nmf_1.Rdata")
save.image(res5,file = "tpm/nmf/nmf_2.Rdata")
###Survival analysis for NMF subgroups
library(survival)
library(survminer)
sur_data <- vector(mode = "list")
for (i in 1:13) {
  sur<-read.table(file = paste("H:/tcga/survival/",tcga_list[i],".survival.tsv",sep = ""),sep = "\t",header = T,check.names = F)
  sur_data[[i]]<-sur[sur$sample %in% colnames(nmf_imput),]
}
sur_data <- do.call("rbind",sur_data)
###
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
###Single cancer type
###Add single cell type group information
nmf_sur_df <- vector(mode = "list")
sur_data <- vector(mode = "list")
for (i in 1:13) {
  sur<-read.table(file = paste("H:/tcga/survival/",tcga_list[i],".survival.tsv",sep = ""),sep = "\t",header = T,check.names = F)
  sur$group <- i
  sur_data[[i]]<-sur[sur$sample %in% colnames(nmf_imput),]
}
sur_data <- do.call("rbind",sur_data)
##single cancer type survival analysis for NMF subgroups
lr_p <- vector()
for (i in 1:13) {
  pdf(file = paste("D:/tcga/model/reconstruct_tpm/tpm/nmf/sur5",i,".pdf",sep = ""),width = 12,height = 10)
  dats<-sur_data[sur_data$group == i,]
  
  fit <- survfit(Surv(OS.time,OS) ~ group2,  
                 data = dats)
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
###
 lr_p<-append(lr_p,NA,after = 3)
###
