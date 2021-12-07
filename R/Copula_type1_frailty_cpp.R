
rm(list=ls())
library(reda)
library(reReg)
library(tidyverse)
library(VineCopula)
require(foreach)
require(pcaPP)
require(doMC)
require(ipred)
require(survival)
require(copula)
require(prodlim)
require(doParallel)
require(Rcpp)
require(RcppArmadillo)
require(frailtypack)
#sourceCpp("src/copula_joint_model.cpp")
#sourceCpp("src/copula_dep.cpp")
#sourceCpp("src/frailty_reg.cpp")
source("R/auc.curve.R")
require(tdROC)
require(modeest)
source("R/hpd.R")


####################################################################
#####the function is for log likelihood function of the data########
### function argurments: #######
##b1: the parameter for x1 in the hazard funtion for the terminal event;
##b2: the parameter for x2 in the hazard funtion for the terminal event;
##b3: the parameter for x1 in the hazard funtion for the recurrent event;
##b4: the parameter for x2 in the hazard funtion for the recurrent event;
##b01, b02: baseline hazards for the terminal event and the recurrent events, respectively;
##phi: random effect term denoted as w in paper; {W_i}
##theta: correlation parameter for copula function; {gamma_ij}
## Mu: theta_mu average level correlation
##x1, x2: two covariates in the hazard function for the terminal event and the recurrent event;
##y1, y2: observed time for the terminal event and the recurrent event;
##s1, s2: failing indicator for the terminal event and the recurrent event;
##id: subject id;
##N: the sample size of subjects;
####################################################################
####################################################################

######################################################
######   generate data and fit the model   ##########
######################################################
simu_data<-function(iseed,t,n,p1,p2,tau,J,gamma){
        set.seed(iseed)
        data_example<-list()
        t<-t
        n<-n ## unique number of ID
        p1<-p1 ## treatment assignment
        x1<-rnorm(n,1,1)
        #x2<-rnorm(n,0,2)
        tau<-tau
        J<-J ## upper for recurrent event number each person
        omega<-rgamma(1,1/tau,1/tau)
        group<-as.numeric(rbinom(n,1,p1))
        gamma<-rep(gamma,n)
        lambda0<-0.9  ## assuming one piece constant
        Beta<-matrix(rnorm(2000,0,1),nrow=1000,ncol=2)
        C<-round(as.numeric(90*rexp(n,lambda0*omega*exp(mean(Beta[,1])*x1))))+1
        C<-ifelse(C>=365,365,C)
        #C<-round(as.numeric(100*rexp(n,exp((x1+x2)))))+1
        #C<-round(as.numeric(rnorm(n,200,10)))
        #C<-rep(365,n)
        id<-seq(1,n,by=1)


        for (i in 1:n){
                j<-1
                r_new<-0
                #gap_time<-round(100*rexp(1,lambda0*omega*exp(group[i]*gamma[i])))+1
                gap_time<-round(60*rexp(1,lambda0*omega*exp(group[i]*gamma[i]
                                                            +mean(Beta[,1])*x1)))+1
                s2<-ifelse(gap_time>=C[i]|gap_time>=t,0,1)
                repeat{
                        if (j>=2)
                                #{r_new<-round(100*rexp(1,lambda0*omega*exp(group[i]*gamma[i])))+1}
                        {r_new<-round(100*rexp(1,lambda0*omega*exp(group[i]*gamma[i]
                                                                   +mean(Beta[,1])*x1)))+1}
                        if (((sum(gap_time)+r_new)<t | (sum(gap_time)+r_new)<C[i]) & length(gap_time)<J) {gap_time=c(gap_time,r_new)
                        s2=c(s2,1)}
                        else if (((sum(gap_time)+r_new)>C[i] |(sum(gap_time)+r_new)>t)  & length(gap_time)<=J+1){
                                gap_time<-c(gap_time,abs(min(t,C[i])-sum(gap_time)))
                                #gap_time=gap_time[-1]
                                #s2<-s2[-1]
                                s2<-c(s2,0)
                                break
                        }
                        else if (((sum(gap_time)+r_new)>C[i]|(sum(gap_time)+r_new)>t) & (sum(gap_time)<min(C[i],t))){
                                gap_time<-c(gap_time,min(t,C[i])-sum(gap_time))
                                #gap_time=gap_time[-1]
                                #s2<-s2[-1]
                                s2<-c(s2,0)
                                break
                        }
                        else if(length(gap_time)>=J+1){
                                gap_time=gap_time[-1]
                                s2=s2[-1]
                                break
                        }
                        j=j+1
                }




                data_example[[i]]<-cbind(id=id[i],gap_time=gap_time/365,group=group[i],lambda0=lambda0,c=C[i]/365,
                                         gamma=gamma[i],event=s2,t=t,tau=tau,days=cumsum(gap_time)/365,x1=x1[i],y1=max(cumsum(gap_time))/365,
                                         ni=1:length(gap_time),maxni=length(gap_time),s1=c(rep(0,length(gap_time)-1),rbinom(1,1,p2)))
        }
        data_example<-data.frame(do.call("rbind",data_example))
        data_example<-data_example[data_example$gap_time!=0,]
        data_example$days<-ifelse(data_example$days>data_example$c,
                                  data_example$c,data_example$days)

        data_example$event<-ifelse(data_example$days<data_example$c|data_example$days<data_example$t,
                                   1,0)
        data_example$event<-ifelse(data_example$days>=data_example$c|data_example$days>=data_example$t,
                                   0,1)
        data_example$s2<-data_example$event
        data_example$y2<-data_example$gap_time
        return(data_example)
}
simulate.data<-function(idx,nsim,n.sim,kk,b1,b2,b3,b4){
        n_patient_s=rep(rep(c(40,80,100),each=9),3) ##possible values for sample size, which is set as 100 ,200 or 400##
        theta_s=rep(c(1,1,1,1.5,1.5,1.5,2,2,2),9)              ##possible values for theta, correlation parameter ##########
        theta=theta_s[idx]              ##for scenario idx=4, theta is 1.5 ##########
        n_patient=n_patient_s[idx]      ##for scenario idx=4, sample size is 100###
        J_s=rep(c(1000,1000,1000),27)    ##for scenario idx=4, J is 100###
        J=J_s[idx]
        UU_s=rep(c(1,0.6,0.4),27)  ## failure rate
        UU=UU_s[idx]
        rho_s=rep(c(0.1,0.3,0.5),each=27)
        rho=rho_s[idx]
        k_s=rep(c(2,4,6),27) ##possible values for number of types of recurrent events ##########
        k=k_s[idx]
        set.seed(n.sim*100)  ###set the seed for simulation ###
        mu=sample(theta_s,k)
        print(mu)
        id=rep(1:n_patient)
        sd1=rnorm(k,0.5,0.001)
        sd2=rnorm(k,0.5,0.001)
        phi<-x1<-x2<-matrix(rep(NA,n_patient*k),nrow = n_patient,ncol = k)
        phi[,kk]=rnorm(n_patient,0,sd1[kk])
        phi[,kk]=scale(phi[,kk],scale=F)  ###random effect w in the paper##
        x1[,kk]=rbinom(n_patient,1,sd1[kk])###a binary covariate########
        x2[,kk]=rnorm(n_patient,0,sd2[kk]) ###a cont covariate############
        b01=rnorm(k,1,sd1[kk])
        b02=rnorm(k,1,sd1[kk])
        #b1=rnorm(k,b1,sd1[kk])
        b1=b1
        b2=rnorm(k,b2,sd2[kk])
        #b3=rnorm(k,b3,sd1[kk])
        b3=b3
        b4=rnorm(k,b4,sd2[kk])
        sigma_epsilon<-0.01
        lam1=exp(b01[kk]+phi[id,kk]+x1[id,kk]*b1[kk]+x2[id,kk]*b2[kk]) ###hazard for the terminal event###
        lam2=exp(b02[kk]+phi[id,kk]+x1[id,kk]*b3[kk]+x2[id,kk]*b4[kk])  ###hazard for the recurrent event###
        u=runif(n_patient)
        t1=-log(u)/lam1    ###generate terminal event###
        c1=runif(n_patient,0,UU)
        y1=pmin(t1,c1,0.5) ###observed time####
        s1=(t1<=c1) ###failing indicator####
        theta_predict=matrix()
        theta=vector()
        main=list()

        gamma=matrix(NA,nrow=n_patient, ncol=J+1)
        for(n_i in 1:n_patient){
                t2=0
                j=1
                gamma[n_i,1]=rnorm(1,0,sqrt(sigma_epsilon/(1-rho^2)))
                theta[kk]=mu[kk]*exp(gamma[n_i,1])
                repeat{
                        if(j>=2) gamma[n_i,j]=rho*gamma[n_i,j-1]+rnorm(1,0,sqrt(sigma_epsilon))

                        theta[kk]=mu[kk]*exp(gamma[n_i,j])

                        if(y1[n_i]==t1[n_i])
                        {
                                w=runif(1)
                                v=((w^(-theta[kk]/(theta[kk]+1))-1)*(u[n_i]^(-theta[kk]))+1)^(-1/theta[kk])
                                #v=-1/theta*log(1+w*(1-exp(-theta))/(w*(exp(-theta*u[n_i])-1)-exp(-theta*u[n_i])))
                                r_new=-log(v)/(lam2[n_i])
                        }

                        else if(s1[n_i]==0){
                                w=exp(-lam1[n_i]*y1[n_i])
                                uv=rCopula(J,claytonCopula(theta[kk]))
                                v_sample=which(uv[,1]<=w)
                                if (length(v_sample)==0 )  uv=rCopula(J,claytonCopula(theta[kk]))

                                v=sample(uv[uv[,1]<=w,2],1,replace=F)
                                r_new=-log(v)/(lam2[n_i])
                                # w=runif(1)
                                # v=((w^(-theta/(theta+1))-1)*(u[n_i]^(-theta))+1)^(-1/theta)
                                # #v=-1/theta*log(1+w*(1-exp(-theta))/(w*(exp(-theta*u[n_i])-1)-exp(-theta*u[n_i])))
                                # r_new=-log(v)/(lam2[n_i])
                        }
                        if(sum(t2)+r_new<=y1[n_i]&length(t2)<J+1) t2=c(t2,r_new)
                        else if (sum(t2)+r_new>y1[n_i]&length(t2)<J+1) {
                                t2=c(t2,y1[n_i]-sum(t2))
                                t2=t2[-1]
                                s2=c(rep(1,length(t2)-1),0)
                                break
                        }
                        else if(length(t2)>=J+1){
                                t2=t2[-1]
                                s2=c(rep(1,length(t2)-1),1)
                                break
                        }
                        j=j+1

                }
                s3=1-s2
                s3[length(s2)]=(s1[n_i])
                main[[n_i]]=cbind(id=n_i,y2=t2,s1=s1[n_i],s2,s3,group=x1[n_i,kk],x2=x2[n_i,kk],y1=y1[n_i],
                                  ni=1:length(t2),maxni=length(t2),theta=theta[kk])
        }

        main=data.frame(do.call("rbind",main))


        return(main)
}
################# estimation in JHCM #############
para_est_JHCM<-function(b01.start,b02.start,b1.start,b2.start,b3.start,b4.start,sigma.start,theta.start,
                        n_sample,n_burn,thin,n.sim,data){
        phi=rnorm(length(unique(data$id)),0,0.5)
        phi=scale(phi,scale=F)  ###random effect w in the paper##
        slice=thin*((n_burn/thin):(n_sample/thin))
        b1_sample=rep(b1.start,n_sample)
        b2_sample=rep(b2.start,n_sample)
        b3_sample=rep(b3.start,n_sample)
        b4_sample=rep(b4.start,n_sample)
        b01_sample=rep(b01.start,n_sample)
        b02_sample=rep(b02.start,n_sample)
        sigma_sample=rep(sigma.start,n_sample)
        theta_sample=rep(theta.start,n_sample)
        scale=c(0.1,0.1,0.1,0.1,0.1,0.1,0.01,0.1) ###step size of MCMC###

        phiN_est_list=MH_joint_model(b1=as.double(b1_sample), b2=as.double(b2_sample),b3=as.double(b3_sample),b4=as.double(b4_sample),
                                     b01=as.double(b01_sample), b02=as.double(b02_sample),
                                     phi=as.double(phi),sigma=as.double(sigma_sample),
                                     theta=as.double(theta_sample), x1=as.double(data$group), x2=as.double(data$x2),
                                     y1=as.double(data$y1),y2=as.double(data$y2), s1=as.integer(data$s1),
                                     s2=as.integer(data$s2), id=as.integer(data$id-1), no_pat_p=as.integer(length(unique(data$id))),
                                     n_sample=as.integer(n_sample), N_p=as.integer(length(data$y1)),scale=as.double(scale))

        param.trace<-cbind(phiN_est_list$b01[slice],phiN_est_list$b1[slice],phiN_est_list$b2[slice],phiN_est_list$b02[slice],
                           phiN_est_list$b3[slice],phiN_est_list$b4[slice],phiN_est_list$sigma[slice],phiN_est_list$theta[slice])
        beta_est<-apply(param.trace,2,mean)
        beta_se<-apply(param.trace,2,sd)
        beta_upper<-beta_est+1.96*beta_se
        beta_lower<-beta_est-1.96*beta_se
        est<-as.data.frame(rbind(beta_est,beta_se, beta_upper,beta_lower))
        colnames(est)<-c("b01","b1","b2","b02","b3,","b4","sigma","theta")

        return(list(param.trace=param.trace,estimation=est,MCMC=phiN_est_list))

}
################# estimation in TVJHCM #############

para_est_TVJHCM<-function(b01.start,b02.start,b1.start,b2.start,b3.start,b4.start,sigma.start,theta.start,
                          rho.start,sigma_epsilon.start,gamma.start,
                          n_sample,n_burn,thin,n.sim,data){
        phi=rnorm(length(unique(data$id)),0,0.5)
        phi=scale(phi,scale=F)  ###random effect w in the paper##
        slice=thin*((n_burn/thin):(n_sample/thin))
        b1_sample=rep(b1.start,n_sample)
        b2_sample=rep(b2.start,n_sample)
        b3_sample=rep(b3.start,n_sample)
        b4_sample=rep(b4.start,n_sample)
        b01_sample=rep(b01.start,n_sample)
        b02_sample=rep(b02.start,n_sample)
        sigma_sample=rep(sigma.start,n_sample)
        rho_sample=rep(rho.start,n_sample)
        theta_sample=rep(na.omit(c(t(gamma.start))),n_sample)
        sigma_epsilon_sample=rep(sigma_epsilon.start,n_sample)
        mu_sample=rep(theta.start,n_sample)

        scale=c(0.1,0.1,0.1,0.1,0.1,0.1,0.01,0.1,0.01)

        result=MH_dep(b1=as.double(b1_sample), b2=as.double(b2_sample),b3=as.double(b3_sample),b4=as.double(b4_sample),
                      b01=as.double(b01_sample), b02=as.double(b02_sample),
                      phi=as.double(phi),sigma=as.double(sigma_sample),mu=as.double(mu_sample),
                      theta=as.double(theta_sample),
                      rho=as.double(rho_sample), sigma_epsilon=as.double(sigma_epsilon_sample), x1=as.double(data$group), x2=as.double(data$x2),
                      y1=as.double(data$y1),y2=as.double(data$y2), s1=as.integer(data$s1),
                      s2=as.integer(data$s2), id=as.integer(data$id-1), no_pat_p=as.integer(length(unique(data$id))),
                      n_sample=as.integer(n_sample),  N_p=as.integer(length(data$y1)),scale=as.double(scale),ni=as.integer(data$ni),maxni=as.integer(data$maxni))

        b1_est=mean(result$b1[slice])
        b2_est=mean(result$b2[slice])
        b3_est=mean(result$b3[slice])
        b4_est=mean(result$b4[slice])
        b01_est=mean(result$b01[slice])
        b02_est=mean(result$b02[slice])
        rho_est=mean(result$rho[slice])
        sigma_epsilon_est=mean(result$sigma_epsilon[slice])
        sigma_est=mean(result$sigma[slice])
        mu_est=mean(result$mu[slice])

        param.trace<-cbind(result$b01[slice],result$b1[slice],result$b2[slice],result$b02[slice],
                           result$b3[slice],result$b4[slice],result$sigma[slice],result$sigma_epsilon[slice],result$mu[slice]
                           ,result$rho[slice])
        beta_est<-apply(param.trace,2,mean)
        beta_se<-apply(param.trace,2,sd)
        beta_upper<-beta_est+1.96*beta_se
        beta_lower<-beta_est-1.96*beta_se
        est<-as.data.frame(rbind(beta_est,beta_se, beta_upper,beta_lower))
        colnames(est)<-c("b01","b1","b2","b02","b3,","b4","sigma","sigma_epsilon","theta_mu","rho")

        return(list(param.trace=param.trace,estimation=est,MCMC=result))

}


############## Aanalysis ###########################################
#################################################################################
######   generate data and fit the model from joint model ##########
#################################################################################

type1<-function(duration,sample_size,p1,p2,tau,J,gamma,sim_num,iseed2){
        output<-list()
        for(iseed in 1:((sim_num)/100)){
                if (iseed>=13){iseed=iseed+1}
                output[[iseed]]<-para_est_JHCM(1,1,0,0,0,0,0.1,1,5000,100,20,iseed,
                                               simu_data(iseed*iseed2,duration/duration,sample_size,p1,p2,tau,J,gamma))$param.trace[,5]
        }

        return(output)

}



result<-list()
tau=1
sample_size=40
### tau=1; gamma=0 ####
result<-type1(365,sample_size,0.5,0.4,tau,30,0,1,1)
#setwd("/Users/mengluliang/Dropbox/biostatisitcs_PSU/2021summerintern/BMS/03results")
saveRDS(result, file = paste0("res_frai_",sample_size,"_",tau,"_cop.RDS"))







