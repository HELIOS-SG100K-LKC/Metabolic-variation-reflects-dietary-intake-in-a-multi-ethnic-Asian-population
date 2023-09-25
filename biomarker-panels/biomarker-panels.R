# To develop dietary biomarker panels and predict intakes
library(tidyverse)
library(stats)
library(glmnet)
library(glue)

# Read files of example datasets (100 participants, 10 foods, 10 metabolites)
data_diet=read.csv("./data/data_diet.csv",row.names=1)
data_met=read.csv("./data/data_met.csv",row.names=1)
data_meta=read.csv("./data/data_meta.csv",row.names=1)

# Merge
diet_met=merge(data_diet,data_met,by="Barcode")
diet_met_meta=merge(data_meta,diet_met,by="Barcode")
rm(diet_met)

# Filter and transform using F1 as example food
data_F1=diet_met_meta %>%
  dplyr::filter(F1<quantile(F1,c(.995),na.rm=TRUE)) %>%
  mutate(F1_log=log1p(F1)) %>%
  dplyr::select(c(age,sex,ethnicity,SSID,group,TYE,incomecat,F1_log,M1:M10)) %>%
  mutate(across(starts_with('M'),~scale(.)))

# Split into train and test sets
train=data_F1 %>%
  filter(group==2) %>%
  drop_na(TYE,incomecat,ethnicity,age) %>%
  select(-c(SSID,group))
test=data_F1 %>%
  filter(group==1) %>%
  drop_na(TYE,incomecat,ethnicity,age) %>%
  select(-c(SSID,group))

train$sex=as.factor(train$sex)
train$ethnicity=as.factor(train$ethnicity)
train$incomecat=as.factor(train$incomecat)

# Convert into matrix
x=model.matrix(train$F1_log~.,train)
y=train$F1_log

# Elastic net penalisation
Penalty_vector=rep(1,ncol(x))
Penalty_vector[1:9]=0

set.seed(888)
# Alpha values scan - level 1
Alpha_values=seq(from=0,to=1,by=0.1)
Min_MSE_values=c()
for (alpha_loop in 1:length(Alpha_values) ) {
  current_alpha_value=0;                current_alpha_value=Alpha_values[alpha_loop]
  cvfit1=cv.glmnet(x,y,penalty.factor=Penalty_vector,family="gaussian",type.measure="dev",standardize=FALSE,alpha=current_alpha_value,nfolds=10,parallel=FALSE)

  Min_MSE_values=append(Min_MSE_values,min(cvfit1$cvm))
}
L1_minIndex=0;           L1_minIndex=which(Min_MSE_values==min(Min_MSE_values) )
L1_min_alpha=0;          L1_min_alpha=Alpha_values[L1_minIndex]
L2_alpha_start=0;        L2_alpha_start=L1_min_alpha[1]-0.05
L2_alpha_start=max(L2_alpha_start,0)
L2_alpha_stop=0;         L2_alpha_stop = L1_min_alpha[1]+0.05
L2_alpha_stop=min(L2_alpha_stop,1)
# Alpha values scan - level 2, resolution: 0.01
L2_Alpha_values=seq(from=L2_alpha_start,to=L2_alpha_stop,by=0.01)
L2_Min_MSE_values=c()
for (alpha_loop in 1:length(L2_Alpha_values) ) {
  current_alpha_value=0;                current_alpha_value=L2_Alpha_values[alpha_loop]
  cvfit2=cv.glmnet(x,y,penalty.factor=Penalty_vector,family="gaussian",type.measure="dev",standardize=FALSE,alpha=current_alpha_value,nfolds=10,parallel=FALSE)
  L2_Min_MSE_values=append(L2_Min_MSE_values,min(cvfit2$cvm))
}
alpha_train=c(Alpha_values,Min_MSE_values,L2_Alpha_values,L2_Min_MSE_values)
plot(cvfit2)
abline(h=min(cvfit2$cvm)+cvfit2$cvsd[which.min(cvfit2$cvm)],col="blue")
plot(Alpha_values,Min_MSE_values, type="o")
plot(L2_Alpha_values,L2_Min_MSE_values, type="o")
minIndex=0;          minIndex=(L2_Min_MSE_values==min(L2_Min_MSE_values))
chosen_alpha_value=L2_Alpha_values[minIndex]
final_cvfit2_1se=glmnet(x,y,penalty.factor=Penalty_vector,standardize=FALSE,alpha=chosen_alpha_value,lambda=cvfit2$lambda.1se,family="gaussian",parallel=FALSE)
coef_elnet_lambda1se=as.matrix(coef(final_cvfit2_1se,s="lambda.1se"))
colnames(coef_elnet_lambda1se)[1]="coef_elnet_lambda1se"

# Ridge regression
set.seed(888)
cvfit3=cv.glmnet(x,y,penalty.factor=Penalty_vector,family="gaussian",type.measure="dev",standardize=FALSE,alpha=current_alpha_value,nfolds=10,parallel=FALSE)
plot(cvfit3)
abline(h=min(cvfit3$cvm)+cvfit3$cvsd[which.min(cvfit3$cvm)],col="blue")
final_cvfit3_1se=glmnet(x,y,penalty.factor=Penalty_vector,standardize=FALSE,alpha=chosen_alpha_value,lambda=cvfit3$lambda.1se,family="gaussian",parallel=FALSE)
coef_ridge_lambda1se=as.matrix(coef(final_cvfit3_1se,s="lambda.1se"))
colnames(coef_ridge_lambda1se)[1]="coef_ridge_lambda1se"

# Merge
coef=data.frame(coef_elnet_lambda1se,coef_ridge_lambda1se)
coef=coef[-c(1:10),]

# Filter
coef$CHEM_ID=rownames(coef)
coef_fil=coef %>% filter(!coef_elnet_lambda1se<=0) %>%
  arrange(desc(coef_elnet_lambda1se)) %>%
  relocate(CHEM_ID,.before=coef_elnet_lambda1se)

# Derive metabolite panel
lm.formula=as.formula(paste0('F1_log ~ ' , 
                                 paste0(glue("{coef_fil$CHEM_ID}"),collapse="+"),'+age+as.factor(sex)+as.factor(ethnicity)+as.factor(incomecat)+TYE'))
model=lm(lm.formula,data=train)
summary(model)
model1=lm(F1_log~age+as.factor(sex)+as.factor(ethnicity),data=train)
summary(model1)
model2 <- lm(F1_log~age+as.factor(sex)+as.factor(ethnicity)+as.factor(incomecat)+TYE, data=train)
summary(model2)

# Partial/incremental F test
fit=anova(model1,model2)
fit$'Pr(>F)'[2] 
fit=anova(model2,model)
fit$'Pr(>F)'[2] 

# Predict
predict_g=predict(model,test)
d=data.frame(test,predict_g)

# Correlation
cor(d$F1_log,d$predict_g,method=c("spearman"),use="pairwise.complete.obs")

# Derive composite score as a weighted sum (ridge regression-derived coefficient) of metabolites from dietary biomarker panels
temp_condition=data_F1 %>%
  dplyr::select(c(which(colnames(data_F1) %in% coef_fil$CHEM_ID))) %>%
  mutate(across(where(is.numeric),~.x*coef_fil[coef_fil$CHEM_ID==cur_column(),3]))
data_F1$predict_score=rowSums(temp_condition)


