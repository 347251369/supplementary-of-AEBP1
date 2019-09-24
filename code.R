##  2019-02-16 ##
library(survival)
library(glmnet)
library(SIS)
library(risksetROC)

load("J:/R_project/gene_project/GBM9.RData")## load data ##
n= 188; pn=10992; 
## genedata - all the gene data with name ##
gene_name=genedata[[1]]        ## the name of gene ##
gene_total_data=rawdata[,-ind1]   ##   rawdata - the origin data  ##




###  strategy - Lasso ####
###   use total data to do Lasso  ####
cv.fit<-cv.glmnet(t(as.matrix(gene_total_data)),Surv(survTime1,cenor1),family="cox")
lambda_best<-cv.fit$lambda.min
fit<-glmnet(t(as.matrix(gene_total_data)),Surv(survTime1,cenor1),family="cox",lambda=lambda_best)
vars.lasso<-nonzeroCoef(fit$beta)               ### the ID of Lasso gene ###
vars.lasso.name<-gene_name[vars.lasso]
gene_lasso_data=gene_total_data[c(vars.lasso),]




###  strategy - Sis ####
### do sis ###
x=t(gene_total_data)
y=Surv(survTime1,cenor1)
model1<-SIS(x,y,family="cox",penalty="lasso",tune="bic",varISIS="van",seed=41)
vars.sis<-model1$ix
vars.sis.name<-gene_name[vars.sis]
gene_sis_data=gene_total_data[c(vars.sis),]




### strategy - SisLasso ####
###  use the variable of Lasso to do condition Sis ###
###  condition sis ###
index_csis<-(1:length(gene_name))[-vars.lasso] 
pvalue_csis<-rep(0, length(index_csis))
for(i in 1:length(index_csis))
{
  j=index_csis[i]
  gene_csis_data=gene_total_data[c(j,vars.lasso),]
  data1<-list(time=survTime1,status=cenor1,ttt=t(as.matrix(gene_csis_data)))
  singlesur.fit=coxph(Surv(time,status)~ttt,data1)
  coef<-singlesur.fit$coefficients[1]
  se<-sqrt(singlesur.fit$var[1,1])
  z<-abs(coef)/se
  p<-pnorm(z)
  pvalue_csis[i]<-2*(1-p)
}
choose<-ifelse(pvalue_csis<=1/pn,1,0)
new0<-index_csis[subset(1:length(index_csis),choose==1)]
new<-c(vars.lasso,new0)
gene_csis_data<-gene_total_data[new,]         ### the condition sis gene data ###

     ### use condition sis data do Lasso again ###
cv.fit<-cv.glmnet(t(as.matrix(gene_csis_data)),Surv(survTime1,cenor1),family="cox")
lambda_best<-cv.fit$lambda.min
fit<-glmnet(t(as.matrix(gene_csis_data)),Surv(survTime1,cenor1),family="cox",lambda =lambda_best)
vars.sislasso<-new[nonzeroCoef(fit$beta)]           ### the ID of SisLasso gene ###
vars.sislasso.name<-gene_name[vars.sislasso]
gene_sislasso_data=gene_total_data[c(vars.sislasso),]





### strategy - improved Sislasso ####
   ###  add the interaction items of sisLasso ### 
row1=nrow(gene_sislasso_data)
col1=ncol(gene_sislasso_data)
gene_improved_sislasso_data<-matrix(0,nrow=row1+row1*(row1-1)/2,ncol=col1)
for(i in 1:row1)
{
  for(j in 1:col1)
  {
    gene_improved_sislasso_data[i,j]<-gene_sislasso_data[i,j]
  }
}
idx=row1;
for(i in 1:(row1-1))
{
  for(j in (i+1):row1)
  {
    idx=idx+1
    for(k in 1:col1)
    {
      gene_improved_sislasso_data[idx,k]<-gene_sislasso_data[i,k]*gene_sislasso_data[j,k]
    }
  }
}
   ###  get the name of sisLasso and add the interaction name ### 
   ##  names(singlesur.fit.sislasso$coefficients)<-gene_sislasso_name   名字放进去## 
gene_sislasso_name<-gene_name[vars.sislasso]
gene_improved_sislasso_name<-vector(length=row1+row1*(row1-1)/2)
gene_improved_sislasso_name<-as.character(gene_sislasso_name)
idx=row1;
for(i in 1:(row1-1))
{
  for(j in (i+1):row1)
  {
    idx=idx+1
    gene_improved_sislasso_name[idx]<-c(paste(paste(gene_sislasso_name[i],"_", sep=""),gene_sislasso_name[j], sep=""))
  }
}



### execute survival ###
sur.lasso_data<-list(time=survTime1,status=cenor1,ttt=t(as.matrix(gene_lasso_data)))
singlesur.fit.lasso=coxph(Surv(time,status) ~ ttt,sur.lasso_data)
names(singlesur.fit.lasso$coefficients)<-vars.lasso.name
summary(singlesur.fit.lasso)

sur.sis_data<-list(time=survTime1,status=cenor1,ttt=t(as.matrix(gene_sis_data)))
singlesur.fit.sis=coxph(Surv(time,status) ~ ttt,sur.sis_data)
names(singlesur.fit.sis$coefficients)<-vars.sis.name
summary(singlesur.fit.sis)

sur.sislasso_data<-list(time=survTime1,status=cenor1,ttt=t(as.matrix(gene_sislasso_data)))
singlesur.fit.sislasso=coxph(Surv(time,status) ~ ttt,sur.sislasso_data)
names(singlesur.fit.sislasso$coefficients)<-vars.sislasso.name
summary(singlesur.fit.sislasso)

sur.improved_sislasso_data<-list(time=survTime1,status=cenor1,ttt=t(as.matrix(gene_improved_sislasso_data)))
singlesur.fit.improved_sislasso=coxph(Surv(time,status) ~ ttt,sur.improved_sislasso_data)
names(singlesur.fit.improved_sislasso$coefficients)<-gene_improved_sislasso_name
summary(singlesur.fit.improved_sislasso)




#### pretreatment of AUC #####
Q1=quantile(utimes,0.25)
Q3=quantile(utimes,0.75)

eta<-singlesur.fit.lasso$linear.predictor
model.score<-eta
utimes<-unique(survTime1[cenor1==1])
utimes<-utimes[order(utimes)]
## find AUC at unique failure times
AUC.lasso<-rep( NA,length(utimes))
for(j in 1:length(utimes))
{
  out<-CoxWeights(eta,survTime1,cenor1,utimes[j])
  AUC.lasso[j]<-out$AUC
}

eta<-singlesur.fit.sis$linear.predictor
model.score<-eta
utimes<-unique( survTime1[cenor1==1])
utimes<-utimes[ order(utimes)]
## find AUC at unique failure times
AUC.sis<-rep( NA,length(utimes))
for(j in 1:length(utimes))
{
  out<-CoxWeights(eta,survTime1,cenor1,utimes[j])
  AUC.sis[j]<-out$AUC
}
 
eta<-singlesur.fit.sislasso$linear.predictor
model.score<-eta
utimes<-unique(survTime1[cenor1==1])
utimes<-utimes[order(utimes)]
## find AUC at unique failure times
AUC.sislasso<-rep( NA,length(utimes))
for(j in 1:length(utimes))
{
  out<-CoxWeights(eta,survTime1,cenor1,utimes[j])
  AUC.sislasso[j]<-out$AUC
}

eta<-singlesur.fit.improved_sislasso$linear.predictor
model.score<-eta
utimes<-unique(survTime1[cenor1==1])
utimes<-utimes[order(utimes)]
## find AUC at unique failure times
AUC.improved_sislasso<-rep( NA,length(utimes))
for(j in 1:length(utimes))
{
  out<-CoxWeights(eta,survTime1,cenor1,utimes[j])
  AUC.improved_sislasso[j]<-out$AUC
}

#### the integrated AUC to get the concordance index inorder to draw ROC curve ######

eta<-singlesur.fit.lasso$linear.predictor
AUC.lasso1=risksetAUC(Stime=survTime1,status=cenor1, marker=eta, method="Cox", tmax=Q1)
AUC.lasso1$Cindex
AUC.lasso3=risksetAUC(Stime=survTime1,status=cenor1, marker=eta, method="Cox", tmax=Q3)
AUC.lasso3$Cindex


eta<-singlesur.fit.sis$linear.predictor
AUC.sis1=risksetAUC(Stime=survTime1,status=cenor1, marker=eta, method="Cox", tmax=Q1)
AUC.sis1$Cindex
AUC.sis3=risksetAUC(Stime=survTime1,status=cenor1, marker=eta, method="Cox", tmax=Q3)
AUC.sis3$Cindex


eta<-singlesur.fit.sislasso$linear.predictor
AUC.sislasso1=risksetAUC(Stime=survTime1,status=cenor1, marker=eta, method="Cox", tmax=Q1)
AUC.sislasso1$Cindex
AUC.sislasso3=risksetAUC(Stime=survTime1,status=cenor1, marker=eta, method="Cox", tmax=Q3)
AUC.sislasso3$Cindex


eta<-singlesur.fit.improved_sislasso$linear.predictor
AUC.improved_csis_cross1=risksetAUC(Stime=survTime1,status=cenor1, marker=eta, method="Cox", tmax=Q1)
AUC.improved_csis_cross1$Cindex
AUC.improved_csis_cross3=risksetAUC(Stime=survTime1,status=cenor1, marker=eta, method="Cox", tmax=Q3)
AUC.improved_csis_cross3$Cindex


### ROC curve  ####
par(mfrow=c(1,2))    ## 画1行2列的图
eta<-singlesur.fit.lasso$linear.predictor
ROC.CC30=risksetROC(Stime=survTime1,status=cenor1,marker=eta,predict.time=30,method="Cox",xlab="specificity",ylab="sensitivity",main="ROC Curve",lty=1,lwd=2,col=2)

eta<-singlesur.fit.sis$linear.predictor
par(new=TRUE)
ROC.CC30=risksetROC(Stime=survTime1,status=cenor1,marker=eta,predict.time=30,method="Cox",xlab="specificity",ylab="sensitivity",main="ROC Curve",lty=1,lwd=2,col=3)

eta<-singlesur.fit.sislasso$linear.predictor
par(new=TRUE)
ROC.CC30=risksetROC(Stime=survTime1,status=cenor1,marker=eta,predict.time=30,method="Cox",xlab="specificity",ylab="sensitivity",main="ROC Curve",lty=1,lwd=2,col=4)

eta<-singlesur.fit.improved_sislasso$linear.predictor
par(new=TRUE)
ROC.CC30=risksetROC(Stime=survTime1,status=cenor1,marker=eta,predict.time=30,method="Cox",xlab="specificity",ylab="sensitivity",main="ROC Curve",lty=1,lwd=2,col=5)

legend("bottomright",legend=c('CoxLasso','CoxSis','CoxSisLasso','CoxImpSisLasso'),col=c(2,3,4,5),lty=1:1,lwd=2:2,text.width=strwidth("1000000000"),cex=0.68)


#### AUC curve ############
Q1=quantile(utimes,0.25)
Q3=quantile(utimes,0.75)
plot(utimes[utimes<=Q3],AUC.sislasso[utimes<=Q3],type="n",ylim=c(0.65,0.95),xlab="time/day",ylab="AUC",font.lab=2,main="Different Methods with AUC")

lines(utimes[utimes<=Q3], AUC.lasso[utimes<=Q3],lty=1,lwd=2,col=2)
lines(utimes[utimes<=Q3], AUC.sis[utimes<=Q3],lty=1,lwd=2,col=3)
lines(utimes[utimes<=Q3], AUC.sislasso[utimes<=Q3],lty=1,lwd=2,col=4)
lines(utimes[utimes<=Q3], AUC.improved_sislasso[utimes<=Q3],lty=1,lwd=2,col=5)

legend("topright",legend=c('CoxLasso','CoxSis','CoxSisLasso','CoxImpSisLasso'),col=c(2,3,4,5),lty=1:1,lwd=2:2,text.width=strwidth("1000000000"),cex=0.68)