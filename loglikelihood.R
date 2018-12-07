

calc_likelihood=function(W_final,beta_final,Xnew,ynew){
Fnew=Xnew%*%W_final

likelihood_mid=data.frame(y=ynew,p=exp(Fnew%*%beta_final)/(1+exp(Fnew%*%beta_final)))
loglik=likelihood_mid$ynew*log(likelihood_mid$p)+(1-likelihood_mid$ynew)*log(1-likelihood_mid$p)
return(loglik)
}
