source("hiscom_code.R")
source("loglikelihood.R")
source("cv_hiscom.R")

#int0 = 1
#intN = 999

index=c(rep("0",16),rep("1",97))
y=as.numeric(index)
group_info=read.csv(,"mRNA_miRNA_groupinfo.csv",header=T,stringsAsFactor=F)
colnames(group_info)=c("group","variable")
dat_wo_val <- read.csv("data_for_pharaoh_for_pval_below_0.05_relationship.csv",header=T,stringsAsFactor=F)

cv_result=perform_cv(dat_wo_val,y,group_info,lambda_list=seq(0.05,0.1,by=0.005),err_rate=0.01,iter_num=10,cvnum=5)
result=perform_hiscom(dat_wo_val,y,group_info,lambda_B=cv_result$lambda_list[which(cv_result$loglik_list==min(cv_result$loglik_list))[1]],lambda_W=cv_result$lambda,err_rate=0.01,iter_num=10)