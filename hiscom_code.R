###Function ÀÛ¼º
	library(glmnet)

perform_hiscom=function(orig_dat,y,group_info,lambda_B,lambda_W,err_rate=0.0001,iter_num){
	num_samp=nrow(orig_dat)
	num_group=length(unique(group_info$group))
	num_variable=nrow(group_info)
	group_list=unique(group_info$group)
	i=1
	genelist=as.numeric(group_info$group==group_list[i])
	if(length(group_list)>1){
	for(i in 2:length(group_list)){
		genelist[which(group_info$group==group_list[i])]=i
		}
	}


	W_init=matrix(0,nrow=num_variable,ncol=num_group)
	i=1
	a=rnorm(sum(genelist==i))
	if(a[1]<0){a=-a}

	W_init[which(genelist==i),i]=a
	if(length(group_list)>1){
		for(i in 2:length(group_list)){
			a=rnorm(sum(genelist==i))
			if(a[1]<0){a=-a}
			W_init[which(genelist==i),i]=a
		}
	}

	X2=as.matrix(orig_dat)
	X2=scale(X2)
	F_init2=X2%*%W_init
	F_init2=scale(F_init2)
	y2=y
	glmfit=glmnet(F_init2,y2,family="binomial",lambda=lambda_B, alpha=0)
	beta_init=rep(0,ncol(W_init))
	beta=as.numeric(glmfit$beta[,1])

	i=1
	XB_init=X2[,which(genelist==i)]*beta[i]
	if(length(group_list)>1){
		for(i in 2:length(group_list)){
			XB_init=cbind(XB_init,X2[,which(genelist==i)]*beta[i])
		}
	}


	glmfit=glmnet(XB_init,y2,family="binomial",lambda=lambda_W ,alpha=0)



	WW=matrix(0,nrow=num_variable,ncol=num_group)

	i=1
	if(glmfit$beta[which(genelist==i)][1]<0){
	WW[which(genelist==i),i]=-glmfit$beta[which(genelist==i)]
	}else{
	WW[which(genelist==i),i]=glmfit$beta[which(genelist==i)]
	}

	if(num_group>1){
	for(i in 2:num_group){
		if(glmfit$beta[which(genelist==i)][1]<0){
			WW[which(genelist==i),i]=-glmfit$beta[which(genelist==i)]
		}else{
			WW[which(genelist==i),i]=glmfit$beta[which(genelist==i)]
		}
	}
	}
 
err1=mean((W_init-WW)^2)
err2=mean((beta_init-beta)^2)
W_init=WW
beta_init=beta

	for(iter in 1:iter_num){

		F_init2=X2%*%W_init
		F_init2=scale(F_init2)

		glmfit=glmnet(F_init2,y2,family="binomial",lambda=lambda_B,alpha=0)
		beta=as.numeric(glmfit$beta[,1])

		i=1
		XB_init=X2[,which(genelist==i)]*beta[i]
		if(num_group>1){
			for(i in 2:num_group){
			XB_init=cbind(XB_init,X2[,which(genelist==i)]*beta[i])
				}
			}

		glmfit=glmnet(XB_init,y2,family="binomial",lambda=lambda_W,alpha=0)

		WW=matrix(0,nrow=num_variable,ncol=num_group)
	i=1
		if(glmfit$beta[which(genelist==i)][1]<0){
			WW[which(genelist==i),i]=-glmfit$beta[which(genelist==i)]
			}else{
			WW[which(genelist==i),i]=glmfit$beta[which(genelist==i)]
		}

		if(num_group>1){
			for(i in 2:num_group){
				if(glmfit$beta[which(genelist==i)][1]<0){
				WW[which(genelist==i),i]=-glmfit$beta[which(genelist==i)]
				}else{
				WW[which(genelist==i),i]=glmfit$beta[which(genelist==i)]
				}
			}
		}
  
	err1=max(abs(W_init-WW))
	err2=max(abs(beta_init-beta))
	err=max(c(err1,err2))
	if(err1<=err_rate & err2<=err_rate){
	break
	}else{
	W_init=WW
	beta_init=beta
	cat(paste(iter,"iteration :",err,"\n"))

		}
	}
result=list(W=WW,beta=beta,err1=err1,err2=err2)
	return(result)
}

