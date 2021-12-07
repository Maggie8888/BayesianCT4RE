
rm(list=ls())
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
require(reReg)

recurrent.event.plot<-function(idx,n.sim){
n_patient_s=rep(rep(c(100,200,400),each=9),3) ##possible values for sample size, which is set as 100 ,200 or 400##
theta_s=rep(c(1,1,1,1.5,1.5,1.5,2,2,2),9)              ##possible values for theta, correlation parameter ##########
theta=theta_s[idx]              ##for scenario idx=4, theta is 1.5 ##########
n_patient=n_patient_s[idx]      ##for scenario idx=4, sample size is 100###
J_s=rep(c(1000,1000,1000),27)    ##for scenario idx=4, J is 100###
J=J_s[idx]
UU_s=rep(c(1,0.6,0.4),27)  ## failure rate
UU=UU_s[idx]
rho_s=rep(c(0.1,0.3,0.5),each=27)
rho=rho_s[idx]
mu=theta

b01=0.5  ###beta_0D=0.5######
b02=1 ###beta_0R=1########
b1=1             ###beta_D1=1#######
b2=1             ###beta_D2=1#######
b3=2             ###beta_R1=2#######
b4=2             ###beta_R2=2#######
n_sample=1000   ###the MCMC posterior sample size is 1000####
n_burn=200     ###the sample size for burn in period is 200####
thin=10          ###thinning the posterior MCMC sample to lower down the correlation#####
slice=thin*((n_burn/thin):(n_sample/thin))
sigma_epsilon=0.1


set.seed(n.sim*100)  ###set the seed for simulation ###
print(n.sim)
id=rep(1:n_patient)
phi=rnorm(n_patient,0,0.5)
phi=scale(phi,scale=F)  ###random effect w in the paper##
x1=rbinom(n_patient,1,0.5) ###a contiunous covariate########
x2=rnorm(n_patient,0,1) ###a binary covariate############
lam1=exp(b01+phi[id]+x1*b1+x2*b2) ###hazard for the terminal event###
lam2=exp(b02+phi[id]+x1*b3+x2*b4) ###hazard for the recurrent event###
u=runif(n_patient)
t1=-log(u)/lam1    ###generate terminal event###
c1=runif(n_patient,0,UU)
y1=pmin(t1,c1,0.5) ###observed time####
s1=(t1<=c1) ###failing indicator####
theta_predict=matrix()
main=list()

gamma=matrix(NA,nrow=n_patient, ncol=J+1)
for(n_i in 1:n_patient){
        t2=0
        t3=0
        j=1
        gamma[n_i,1]=rnorm(1,0,sqrt(sigma_epsilon/(1-rho^2)))
        theta=mu*exp(gamma[n_i,1])
        repeat{
                if(j>=2) gamma[n_i,j]=rho*gamma[n_i,j-1]+rnorm(1,0,sqrt(sigma_epsilon))

                theta=mu*exp(gamma[n_i,j])

                if(y1[n_i]==t1[n_i])
                {
                        w=runif(1)
                        v=((w^(-theta/(theta+1))-1)*(u[n_i]^(-theta))+1)^(-1/theta)
                        #v=-1/theta*log(1+w*(1-exp(-theta))/(w*(exp(-theta*u[n_i])-1)-exp(-theta*u[n_i])))
                        r_new=-log(v)/(lam2[n_i])
                }

                else if(s1[n_i]==0){
                        w=exp(-lam1[n_i]*y1[n_i])
                        uv=rCopula(J,claytonCopula(theta))
                        v_sample=which(uv[,1]<=w)
                        if (length(v_sample)==0 )  uv=rCopula(J,claytonCopula(theta))

                        v=sample(uv[uv[,1]<=w,2],1,replace=F)
                        r_new=-log(v)/(lam2[n_i])
                        # w=runif(1)
                        # v=((w^(-theta/(theta+1))-1)*(u[n_i]^(-theta))+1)^(-1/theta)
                        # #v=-1/theta*log(1+w*(1-exp(-theta))/(w*(exp(-theta*u[n_i])-1)-exp(-theta*u[n_i])))
                        # r_new=-log(v)/(lam2[n_i])
                }
                if(sum(t2)+r_new<=y1[n_i]&length(t2)<J+1) {t2=c(t2,r_new)
                t3=c(t3,sum(t2))}
                else if (sum(t2)+r_new>y1[n_i]&length(t2)<J+1) {
                        t2=c(t2,y1[n_i]-sum(t2))
                        t2=t2[-1]
                        s2=c(rep(1,length(t2)-1),0)
                        t3=c(t3)
                        break
                }
                else if(length(t2)>=J+1){
                        t2=t2[-1]
                        s2=c(rep(1,length(t2)-1),1)
                        t3=c(t3)
                        break
                }
                j=j+1

        }
        s3=1-s2
        s3[length(s2)]=(s1[n_i])
        main[[n_i]]=cbind(id=n_i,y2=t2,s1=c(rep(0,length(t2)-1),1-s1[n_i]),s2,s3,x1=x1[n_i],x2=x2[n_i],y1=y1[n_i],ni=1:length(t2),maxni=length(t2),t3=t3)
}

main=data.frame(do.call("rbind",main))
main$t3<-ifelse(main$y1<=main$t3,main$y1,main$t3)
reObj <- Recur(time = main$t3, id = main$id, event = main$s2, terminal = main$s1)
return(reObj)
}

p1<-recurrent.event.plot(81,1)

p2<-recurrent.event.plot(1,1)

p3<-recurrent.event.plot(37,1)

par(mfrow=c(1,3))
plot(p2)
plot(p3)
plot(p1)
