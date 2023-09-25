# To assess partial correlation coefficients between diet variables and metabolites controlled for effects of age, sex, ethnicity and batch.
library(psych)
library(corrr)
library(gplots)
library(tidyverse)
library(stringr)

# Read files of example datasets (100 participants, 10 foods, 10 metabolites)
data_diet=read.csv("./data/data_diet.csv",row.names=1)
data_met=read.csv("./data/data_met.csv",row.names=1)
data_meta=read.csv("./data/data_meta.csv",row.names=1)
  
# Merge
diet_met=merge(data_diet,data_met,by="Barcode")
diet_met_meta=merge(data_meta,diet_met,by="Barcode")
rm(diet_met)

# Transform
diet_met_log=diet_met_meta %>%
  select(contains('F'),Barcode,age,sex,ethnicity) %>%
  mutate(across(contains('F'),~log1p(.))) %>%
  rename_all(~stringr::str_replace(.,regex("^F",ignore_case=TRUE),"F_log"))
diet_met_meta=diet_met_meta %>%
  left_join(diet_met_log)
rm(diet_met_log)

# Partial correlation
result_r=partial.r(diet_met_meta,x=c(19:ncol(diet_met_meta)),y=c(2:4),method="spearman")

# Adjust p values
result_p=corr.p(result_r,n=20,adjust="none",ci=FALSE)$p
result_p_fdr=result_p
result_p_fdr[]=p.adjust(result_p,method="BH",n=length(result_p))

# Filter
result_r_fil=as.data.frame(result_r[-c(1:10),-c(11:20)])
result_r_fil=result_r_fil %>%
  select_if(~any(.>=0.15))
result_p_fdr_fil=as.data.frame(result_p_fdr[-c(1:10),-c(11:20)])
result_p_fdr_fil=result_p_fdr_fil %>%
  select_if(~any(.>=0.15))

write.csv(result_r_fil,paste('./results/pcorr_r.csv',sep=''))
write.csv(result_p_fdr_fil,paste('./results/pcorr_fdr_p.csv',sep=''))

# Set sequence for heatmap colour
Breaks=seq(-0.7,0.7,0.01)
colours=bluered(140)

# Heatmap
svg(file="heatmap.svg",width=12, height=17)
heatmap=heatmap.2(as.matrix(result_r_fil[,-c(1)]),scale="none",trace="none",density.info="none",col=colours,margins=c(8,20),
                  keysize=0.5,key.par=list(cex=0.4),
                  Rowv=TRUE,
                  Colv=NULL,
                  lwid=c(1,5),
                  lhei=c(0.7,6),
                  cexRow=0.8, cexCol=0.8,
                  dendrogram="none",breaks=Breaks)
dev.off()

