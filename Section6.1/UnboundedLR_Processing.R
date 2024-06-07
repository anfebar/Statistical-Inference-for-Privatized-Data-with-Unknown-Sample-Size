library(xtable)
### set the working directory to source file location

#file_path = paste("./output/eps1_ep_N1_iteration",1,".RData",sep="")
#load(file_path)
burn_in = 5000
ep_vec = c(.1,1)
beta_index = 2
N_DP_vec = c(0.001,.01,.1,1,10,Inf)
#N_DP_vec = c(Inf)
tables = list()



getSummaries = function(file_path){
  load(file_path)
  beta_mean = as.numeric(rowMeans(chain$beta[,-(2:burn_in)]))
  beta_var = as.numeric(apply(chain$beta[,-(1:burn_in)],1,var))
  n_mean = mean(chain$n[-(1:burn_in)])
  n_var = var(chain$n[-(1:burn_in)])
  tau_mean = mean(chain$tau[-(1:burn_in)])
  tau_var = var(chain$tau[-(1:burn_in)])
  return(list(beta_mean=beta_mean,beta_var=beta_var,n_mean=n_mean,n_var=n_var,
              tau_mean = tau_mean, tau_var=tau_var))
}


for(ep in ep_vec){
  print(c("epsilon: ",ep))
Beta_M = rep(NaN,length(N_DP_vec))
Beta_M_SE = rep(NaN,length(N_DP_vec))
Beta_V = rep(NaN,length(N_DP_vec))
Beta_V_SE = rep(NaN,length(N_DP_vec))      

N_M = rep(NaN,length(N_DP_vec))
N_M_SE = rep(NaN,length(N_DP_vec))
N_V = rep(NaN,length(N_DP_vec))
N_V_SE = rep(NaN,length(N_DP_vec))

Tau_M = rep(NaN,length(N_DP_vec))
Tau_M_SE = rep(NaN,length(N_DP_vec))
Tau_V = rep(NaN,length(N_DP_vec))
Tau_V_SE = rep(NaN,length(N_DP_vec))
                          
             
for(e in 1:length(N_DP_vec)){
  print(e)
  N_DP = N_DP_vec[e]
  
  
  
  file_prefix=""
  if(N_DP==Inf){
    file_prefix = paste("./outputs/eps",ep,"_Bounded_iteration",sep="")
  }else{
    file_prefix = paste("./outputs/eps",ep,"_ep_N",N_DP,"_iteration",sep="")
  }
  
  

  beta_means = matrix(rep(NaN,100*3),nrow=3,ncol=100)
  beta_vars = matrix(rep(NaN,100*3),nrow=3,ncol=100)
  n_means = rep(NaN,100)
  n_vars = rep(NaN,100)
  tau_means = tau_vars = rep(NaN,100)
  
  for(i in 1:100){
    print(i)
    file_path = paste(file_prefix,i,".RData",sep="")
    
    tryCatch({
    summaries = getSummaries(file_path)
      beta_means[,i] = summaries$beta_mean
      beta_vars[,i] = summaries$beta_var
      n_means[i] = summaries$n_mean
      n_vars[i] = summaries$n_var
      tau_means[i] = summaries$tau_mean
      tau_vars[i] = summaries$tau_var
    },error=function(e) {
      message('An Error Occurred')
      print(e)
    },
    #if a warning occurs, tell me the warning
    warning=function(w) {
      message('A Warning Occurred')
      print(w)
      return(NA)
    }
    )
    #print(summaries)
  }
  
  
  Beta_M[e] = mean(beta_means[beta_index,],na.rm=TRUE)
  Beta_M_SE[e] = sqrt(var(beta_means[beta_index,],na.rm=TRUE)/length(n_means))### monte carlo standard errors
  Beta_V[e] = mean(beta_vars[beta_index,],na.rm=TRUE)
  Beta_V_SE[e] = sqrt(var(beta_vars[beta_index,],na.rm=TRUE)/length(n_means))
  N_M[e] = mean(n_means,na.rm=TRUE)
  N_M_SE[e] = sqrt(var(n_means,na.rm=TRUE)/length(n_means))
  N_V[e] = mean(n_vars,na.rm=TRUE)
  N_V_SE[e] = sqrt(var(n_vars,na.rm=TRUE)/length(n_means))
  Tau_M[e] = mean(tau_means,na.rm=TRUE)
  Tau_M_SE[e] = sqrt(var(tau_means,na.rm=TRUE)/length(n_means))
  Tau_V[e] = mean(tau_vars,na.rm=TRUE)
  Tau_V_SE[e] = sqrt(var(tau_vars,na.rm=TRUE)/length(n_means))
}

#plot(log(N_DP_vec),Beta_M,main=paste("ep=",ep))
#plot(log(N_DP_vec),Beta_V,main=paste("ep=",ep))
#plot(log(N_DP_vec),Tau_M,main=paste("ep=",ep))
#plot(log(N_DP_vec),Tau_V,main=paste("ep=",ep))

#plot(log(N_DP_vec),N_V,main=paste("ep=",ep))
#Beta_M_SE
#Beta_V_SE
#Tau_M_SE
#Tau_V_SE
#N_V_SE



df = as.data.frame(round(1000*rbind(Beta_M, Beta_V, Tau_M, Tau_V, N_M, N_V))/1000, stringsAsFactors = FALSE)
names(df) <- N_DP_vec
#df

tables[[paste(ep,"_estimates",sep="")]]=df

df_SE = as.data.frame(round(1000*rbind(Beta_M_SE, Beta_V_SE, Tau_M_SE, Tau_V_SE, N_M_SE, N_V_SE))/1000, stringsAsFactors = FALSE)
names(df_SE) <- N_DP_vec
#df_SE
tables[[paste(ep,"_errors",sep="")]]=df_SE
}
tables

#### for latex output

xtable(tables$'0.1_estimates',digits=3,align="lrrrrrr")

xtable(tables$'0.1_errors',digits=3,align="lrrrrrr")

xtable(tables$'1_estimates',digits=3,align="lrrrrrr")

xtable(tables$'1_errors',digits=3,align="lrrrrrr")

########################################################
ep=1
N_DP = .01
i=1
file_path2 =paste("./outputs/eps",ep,"_ep_N",N_DP,"_iteration",i,".RData",sep="")
load(file_path2)
pdf("./figures/beta_1_01.pdf",width=5,height=5)
matplot(t(chain$beta_chain), type = "l",ylim=c(-3,3),xlab="iteration",ylab="beta chain");
abline(h=0,col="black")
abline(h=-1,col="red")
abline(h=1,col="green")
dev.off()
pdf("./figures/n_1_01.pdf",width=5,height=5)
plot(chain$n_chain,type="l",xlab="iteration",ylab="N chain")
abline(h=1000)
dev.off()
##########