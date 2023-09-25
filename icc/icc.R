# To calculate intraclass correlation (ICC) between visits
library(tidyverse)
library(irr)

# Read files
data_icc=read.csv("./data/data_icc.csv",row.names=1)

# ICC for predicted intakes
data_intake=data_icc %>%
  dplyr::select("Barcode","Visit_number","predict_g") %>%
  arrange(Barcode) %>%
  pivot_wider(names_from="Visit_number",values_from="predict_g") %>%
  mutate(rater=matrix("rater",nrow=195,ncol=1)) %>%
  mutate(num=matrix(1:195,nrow=195,ncol=1)) %>%
  unite("rater_num",rater:num,sep="") %>%
  relocate(rater_num,.before=2) %>%
  select(-c(Barcode,rater_num))
data_intake_icc=icc(data_icc,model,"two-way",type="consistency",unit="single",r0=0,conf.level=0.95)

# ICC for reported intakes
data_intake=data_icc %>%
  dplyr::select("Barcode","Visit_number","orange_g_day") %>%
  arrange(Barcode) %>%
  pivot_wider(names_from="Visit_number",values_from="orange_g_day") %>%
  mutate(rater=matrix("rater",nrow=195,ncol=1)) %>%
  mutate(num=matrix(1:195,nrow=195,ncol=1)) %>%
  unite("rater_num",rater:num,sep="") %>%
  relocate(rater_num,.before=2) %>%
  select(-c(Barcode,rater_num))
data_intake_icc=icc(data_icc,model,"two-way",type="consistency",unit="single",r0=0,conf.level=0.95)

# ICC for metabolite scores
data_score=data_icc %>%
  dplyr::select("Barcode","Visit_number","predict_score") %>%
  arrange(Barcode) %>%
  pivot_wider(names_from="Visit_number",values_from="predict_score") %>%
  mutate(rater=matrix("rater",nrow=195,ncol=1)) %>%
  mutate(num=matrix(1:195,nrow=195,ncol=1)) %>%
  unite("rater_num",rater:num,sep="") %>%
  relocate(rater_num,.before=2) %>%
  select(-c(Barcode,rater_num))
data_intake_icc=icc(data_icc,model,"two-way",type="consistency",unit="single",r0=0,conf.level=0.95)

