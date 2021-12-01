#include <string>
#include <RcppArmadillo.h>
#include <math.h>
#include <Rmath.h>
#include <Rcpp.h>
#include <list> 
#include <iterator> 
//#include <bits/stdc++.h> 
#include <vector> 
#include <random>
#include <algorithm>
#include <iostream>
#include <cstdlib>
#include <list> 
#include <stdlib.h>
#include <type_traits>
using namespace std; 
using namespace arma;
using namespace Rcpp;
#include <iostream>

//using std::swap;
//using std::cout;
//using std::endl;


// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

NumericMatrix zero_mat(int rows, int cols)
{
        int i, j;
        NumericMatrix matrix;
        for (i = 0; i < rows; i++)
        {
                for (j = 0; j < cols; j++) 	
                        matrix(i,j)=0;
        }
        return(matrix);
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

NumericMatrix ident_mat(int dim)
{
        int i, j;
        NumericMatrix matrix;
        for (i = 0; i < dim; i++)
        {
                for (j = 0; j < dim; j++)
                { if (i!=j) 
                        matrix(i,j)=0;
                else matrix(i,j)=1;
                }
        }
        
}


// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

NumericVector center_dep(NumericVector x,int n_region){
        double sum=0.0, avg;
        for (int i=0;i<n_region;i++){
                sum=sum+x[i];
        }
        avg=sum/n_region;
        for (int j=0;j<n_region;j++){
                x[j]=x[j]-avg;
        }
        return x;
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

NumericMatrix transpose (NumericMatrix matrix,int rows, int cols)
{
        int  i, j;
        NumericMatrix x;
        for (i = 0; i < rows; i++) {
                for ( j = 0; j < cols; j++) {
                        x(j,i) = matrix(i,j);
                }
        }
        return(x);
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
NumericMatrix matpr (NumericMatrix a,NumericMatrix b, int m, int p, int n)
{
        int i, j, k;
        double sum;
        NumericMatrix c;
        for (i = 0; i < m ; i++) {
                for(j = 0; j < n; j++) {
                        sum = 0 ;
                        for (k = 0; k < p; k++) {
                                sum += a(i,k) * b(k,j);
                        }
                        c(i,j) = sum;
                }
        }
        return(c);
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

NumericMatrix matsum(NumericMatrix a,NumericMatrix b, int rows, int cols)
{
        int i, j;
        NumericMatrix matrix;
        for (i = 0; i < rows; i++) {
                for (j = 0; j < cols; j++) 
                        matrix(i,j) = a(i,j) + b(i,j);
        }
        return(matrix);
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

double SumOfMatrix(NumericMatrix  mat, int n, int m)
{
        
        double sum = 0;
        int i, j;
        
        for (i = 0; i < n; ++i) {
                for (j = 0; j < m; ++j) {
                        sum += mat(i,j);
                }
        }
        
        return(sum);
        
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

NumericMatrix vectormat (NumericVector vector, int nrow, int ncol)
{
        int i, j;
        NumericMatrix matrix;
        for (i = 0; i < nrow; i++) {
                for (j = 0; j < ncol; j++) 
                        matrix(i,j)= vector[j*nrow+i];
        }
        return(matrix);
}



// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

double LLi(double tau,NumericVector lambda0,NumericVector t, 
           NumericVector z,NumericMatrix x,NumericMatrix beta,NumericVector c,double gamma){
        NumericMatrix temp,mat,reZ;
        double sum=0;
        int betaRows = std::extent<decltype(beta), 0>::value;
        int betaCols = std::extent<decltype(beta), 1>::value;
        int xRows = std::extent<decltype(x), 0>::value;
        int xCols = std::extent<decltype(x), 1>::value;
       // int tbetaCols = std::extent<decltype(transpose(beta)), 1>::value;
        
        for(int i=0;i<z.size();++i){
                 reZ(i,0) = z[i] * gamma;
        }
        
        int reZRows = std::extent<decltype(reZ), 0>::value;
        int reZCols = std::extent<decltype(reZ), 1>::value;
         temp=matpr(x,beta,xRows,xCols,betaCols);
         mat=matsum(reZ,temp,reZRows,reZCols);
         int matRows = std::extent<decltype(mat), 0>::value;
         int matCols = std::extent<decltype(mat), 1>::value;
         double summ=SumOfMatrix(mat,matRows,matCols);
        for (int j=0; j<t.length();j++){
                if (t[j]<=c[j]){
                        double gamma0=lambda0[j]*t[j];
        sum+=sum+log(tau)+log(lambda0[j])+summ+
                lgamma(1/tau+t.length())-lgamma(1/tau)
                        -(1/tau+t.length())*log(1+tau*gamma0*exp(summ));
                }
        }
        return(sum);
}


// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
double LL(NumericVector id,double TAU,NumericVector lambda0,NumericVector T,
          NumericVector Z,NumericMatrix X,NumericMatrix Beta,int N,NumericVector C,
        double gamma,NumericVector uniqueid,int no_p){
        double sum=0, tau;
        int XCols = std::extent<decltype(X), 1>::value;
            for (int l=0;l<no_p;l++){
                NumericVector z,t,c,temp;
                NumericMatrix x;
               for(int i=0;i<N;++i){
                  t=T[id[i]==uniqueid[l]];
                  z=Z[id[i]==uniqueid[l]];
                  temp=X(id[i]==uniqueid[l],_);
                  for(int m=0;m<temp.length();m++){
                        for (int k=0;k<XCols;++k){
                         x(m,k)=X(id[i]==uniqueid[l],k);
                       }
                  }
                  c=C[id[i]==uniqueid[l]];
                  tau=TAU;
                  sum=sum+LLi(tau,lambda0,t,z,x,Beta,c,gamma);
               }
            }
        return(sum);
}


// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

double LLi2(double tau,NumericVector lambda0,NumericVector t, 
           NumericVector z,NumericMatrix x,NumericMatrix beta,NumericVector c,double gamma,
           double alpha){
        NumericMatrix temp,mat,reZ;
        double sum=0;
        int betaRows = std::extent<decltype(beta), 0>::value;
        int betaCols = std::extent<decltype(beta), 1>::value;
        int xRows = std::extent<decltype(x), 0>::value;
        int xCols = std::extent<decltype(x), 1>::value;
        //int tbetaCols = std::extent<decltype(transpose(beta)), 1>::value;
        
        for(int i=0;i<z.size();++i){
                reZ(i,0) = z[i] * gamma;
        }
        
        int reZRows = std::extent<decltype(reZ), 0>::value;
        int reZCols = std::extent<decltype(reZ), 1>::value;
        temp=matpr(x,beta,xRows,xCols,betaCols);
        mat=matsum(reZ,temp,reZRows,reZCols);
        int matRows = std::extent<decltype(mat), 0>::value;
        int matCols = std::extent<decltype(mat), 1>::value;
        double summ=SumOfMatrix(mat,matRows,matCols);
        for (int j=0; j<t.length();j++){
                if (t[j]<=c[j]){
                        double gamma0=lambda0[j]*t[j];
                        sum+=sum+log(tau)+log(lambda0[j])+summ+
                                lgamma(1/tau+t.length())-lgamma(1/tau)
                                -(1/tau+t.length())*log(1+tau*gamma0*exp(summ));
                }
        }
        return(alpha*sum);
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
double LL2(NumericVector id,double TAU,NumericVector lambda0,NumericVector T,
          NumericVector Z,NumericMatrix X,NumericMatrix Beta,int N,NumericVector C,
          double gamma,double alpha,NumericVector uniqueid,int no_p){
    double sum=0, tau;
    int XCols = std::extent<decltype(X), 1>::value;
    for (int l=0;l<no_p;l++){
        NumericVector z,t,c,temp;
        NumericMatrix x;
        for(int i=0;i<N;++i){
            t=T[id[i]==uniqueid[l]];
            z=Z[id[i]==uniqueid[l]];
            temp=X(id[i]==uniqueid[l],_);
            for(int m=0;m<temp.length();m++){
                for (int k=0;k<XCols;++k){
                    x(m,k)=X(id[i]==uniqueid[l],k);
                }
            }
            c=C[id[i]==uniqueid[l]];
            tau=TAU;
            sum=sum+LLi2(tau,lambda0,t,z,x,Beta,c,gamma,alpha);
        }
    }
    return(sum);
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]


double target_omega(NumericVector id,double TAU,NumericVector lambda0,NumericVector T,
                    NumericVector Z,NumericMatrix X,NumericMatrix Beta,int N,NumericVector C,
                    double gamma,int k,NumericVector omega,double x1_new,NumericVector uniqueid,int no_p)
{
        double sum = 0.0, func=0.0;
        int m;
        NumericVector temp=omega;
        for(m=1; m<N; m++)
        {
                if(m==k) {
                        temp[m]=x1_new;
                }
                else
                {
                        temp[m]=omega[m];
                }
        }

        func=func+log(R::dgamma(temp[m],1/TAU,1/TAU,0));
        sum=LL(id,TAU,lambda0,T,Z,X,Beta,N,C,
        gamma,uniqueid,no_p)-0.5*func;
        return (sum);
}


// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]


double target_omega2(NumericVector id,NumericVector id2,double TAU,NumericVector lambda0,NumericVector T,NumericVector T2,
                    NumericVector Z,NumericVector Z2,NumericMatrix X,NumericMatrix X2,
                    NumericMatrix Beta,int N,int N2,NumericVector C,NumericVector C2,
                    double gamma,int k,NumericVector omega,double x1_new,double alpha,NumericVector uniqueid,NumericVector uniqueid2,
                    int no_p,int no_p2)
{
        double sum = 0.0, func=0.0;
        int m;
        NumericVector temp=omega;
        for(m=1; m<no_p; m++)
        {
                if(m==k) {
                        temp[m]=x1_new;
                }
                else
                {
                        temp[m]=omega[m];
                }
        }
        
        func=func+log(R::dgamma(temp[m],1/TAU,1/TAU,0));
        sum=LL2(id,TAU,lambda0,T,Z,X,Beta,N,C,
               gamma,alpha,uniqueid,no_p)+LL(id2,TAU,lambda0,T2,Z2,X2,Beta,N2,C2,gamma,uniqueid2,no_p2)-0.5*func;
        return (sum);
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

double target_tau(NumericVector omega,int N)
{
        double sum=0.0;
        for(int i=0; i<N; i++){
                sum=sum+omega[i]*omega[i];
        }
        return(sum);
}


// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]


double target_alpha(NumericVector id,NumericVector id2,double TAU,NumericVector lambda0,NumericVector T,
                    NumericVector T2,NumericVector Z,NumericVector Z2,NumericMatrix X,NumericMatrix X2,
                    NumericMatrix Beta,int N,int N2,NumericVector C,NumericVector C2,
                    double gamma,int k,NumericVector alpha,double x1_new,NumericVector uniqueid,NumericVector uniqueid2,
                    int no_p,int no_p2)
{
        double sum = 0.0, func=0.0;
        int m;
        NumericVector temp=alpha;
        for(m=1; m<no_p; m++)
        {
                if(m==k) {
                        temp[m]=x1_new;
                }
                else
                {
                        temp[m]=alpha[m];
                }
        }
        
        func=func+log(R::dunif(temp[m],0,1,0));
        sum=LL2(id,TAU,lambda0,T,Z,X,Beta,N,C,
               gamma,temp[m],uniqueid,no_p)
            //+LL(id2,TAU,lambda0,T2,Z2,X2,Beta,N2,C2,gamma,uniqueid2,no_p2)
            +func;
        return (sum);
}





// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

double rand_gen() {
        // return a uniformly distributed random value
        return ( (double)(rand()) + 1. )/( (double)(RAND_MAX) + 1. );
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

double normalRandom() {
        // return a normally distributed random value
        double v1=rand_gen();
        double v2=rand_gen();
        return cos(2*3.14*v2)*sqrt(-2.*log(v1));
}


// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

double metroplis_omega(NumericVector id,double TAU,NumericVector lambda0,NumericVector T,
                       NumericVector Z,NumericMatrix X,NumericMatrix Beta,int N,NumericVector C,
                       double gamma,int k,NumericVector omega_curr,double x1_prev, double scale,
                       NumericVector uniqueid,int no_p)
{
        double lgratio=0,temp1=0,temp2=0;
        double x1_curr;
        double u=R::runif(0,1);
        double sd = 1.0;
        double mean = 0.0;
        double rnorm = normalRandom()*sd+mean;
        x1_curr=x1_prev+scale*rnorm;
        temp1=target_omega(id,TAU,lambda0,T,Z,X,Beta,N,C,
                           gamma,k,omega_curr,x1_curr,uniqueid,no_p);
        temp2=target_omega(id,TAU,lambda0,T,Z,X,Beta,N,C,
                           gamma,k,omega_curr,x1_prev,uniqueid,no_p);
        lgratio=temp1-temp2;
        
        if (lgratio>=log(u))
        {
                omega_curr[k]=x1_curr;
        }
        else
        {
                omega_curr[k]=x1_prev;
        }
        
        return  omega_curr[k];
}


// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

double metroplis_omega2(NumericVector id,NumericVector id2,double TAU,NumericVector lambda0,NumericVector T,
                        NumericVector T2,NumericVector Z,NumericVector Z2,NumericMatrix X,NumericMatrix X2,
                        NumericMatrix Beta,int N,int N2,NumericVector C,NumericVector C2,
                       double gamma,int k,NumericVector omega_curr,double x1_prev,double alpha, double scale,
                       NumericVector uniqueid,NumericVector uniqueid2,int no_p,int no_p2)
{
        double lgratio=0,temp1=0,temp2=0;
        double x1_curr;
        double u=R::runif(0,1);
        double sd = 1.0;
        double mean = 0.0;
        double rnorm = normalRandom()*sd+mean;
        x1_curr=x1_prev+scale*rnorm;
        temp1=target_omega2(id,id2,TAU,lambda0,T,T2,Z,Z2,X,X2,Beta,N,N2,
                            C,C2,gamma,k,omega_curr,x1_curr,alpha,uniqueid,uniqueid2,no_p,no_p2);
        temp2=target_omega2(id,id2,TAU,lambda0,T,T2,Z,Z2,X,X2,Beta,N,N2,
                            C,C2,gamma,k,omega_curr,x1_prev,alpha,uniqueid,uniqueid2,no_p,no_p2);
        lgratio=temp1-temp2;
        
        if (lgratio>=log(u))
        {
                omega_curr[k]=x1_curr;
        }
        else
        {
                omega_curr[k]=x1_prev;
        }
        
        return  omega_curr[k];
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

double metroplis_alpha(NumericVector id,NumericVector id2,double TAU,NumericVector lambda0,NumericVector T,
                       NumericVector T2,NumericVector Z,NumericVector Z2,NumericMatrix X,NumericMatrix X2,
                       NumericMatrix Beta,int N,int N2,NumericVector C,NumericVector C2,
                       double gamma,int k,NumericVector alpha_curr,double x1_prev,double scale,
                       NumericVector uniqueid,NumericVector uniqueid2,int no_p,int no_p2)
{
        double lgratio=0,temp1=0,temp2=0;
        double x1_curr;
        double u=R::runif(0,1);
        double sd = 1.0;
        double mean = 0.0;
        double rnorm = normalRandom()*sd+mean;
        x1_curr=x1_prev+scale*rnorm;
        temp1=target_alpha(id,id2,TAU,lambda0,T,T2,Z,Z2,X,X2,Beta,N,N2,
                           C,C2,gamma,k,alpha_curr,x1_curr,uniqueid,uniqueid2,no_p,no_p2);
        temp2=target_alpha(id,id2,TAU,lambda0,T,T2,Z,Z2,X,X2,Beta,N,N2,
                           C,C2,gamma,k,alpha_curr,x1_prev,uniqueid,uniqueid2,no_p,no_p2);
        lgratio=temp1-temp2;
        
        if (lgratio>=log(u))
        {
                alpha_curr[k]=x1_curr;
        }
        else
        {
                alpha_curr[k]=x1_prev;
        }
        
        return  alpha_curr[k];
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]


double target_beta(NumericVector id,double TAU,NumericVector lambda0,NumericVector T,
                    NumericVector Z,NumericMatrix X,NumericMatrix Beta,int N,NumericVector C,
                    double gamma,int k,NumericMatrix omega,double n_sample,NumericVector x1_new,
                    NumericVector uniqueid,int no_p)
{
        double sum = 0.0, func=0.0;
        int m;
        NumericMatrix tempp;
        NumericMatrix temp=omega;
        for(m=1; m<n_sample; m++)
        {
                if(m==k) {
                        temp(m,_)=x1_new;
                }
                else
                {
                        temp(m,_)=omega(m,_);
                }
        }
        
        //func=func+log(R::dgamma(temp(m,0),1,1,0));
        for (int s=0; s<temp(m,_).size(); s++){
            tempp(s,0)=temp(m,s);
        }
        
        sum=LL(id,TAU,lambda0,T,Z,X,tempp,N,C,
               gamma,uniqueid,no_p);
        return (sum);
}


// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

NumericVector metroplis_beta(NumericVector id,double TAU,NumericVector lambda0,NumericVector T,
                      NumericVector Z,NumericMatrix X,NumericMatrix Beta,
                      int N,NumericVector C,
                      double gamma,int k,NumericMatrix Beta_curr,NumericVector x1_prev,
                      double n_sample,double scale,NumericVector uniqueid,int no_p)
{
        double lgratio=0, u ,temp1=0, temp2=0;
        NumericVector x1_curr;
        u=R::runif(0,1);
        for(int j=0; j<x1_prev.length(); j++){
              x1_curr[j]=x1_prev[j]+scale*R::rnorm(0,1);  
        }
        
        temp1=target_beta(id,TAU,lambda0,T,Z,X,Beta,N,C,gamma,k,Beta_curr,n_sample,x1_curr,uniqueid,no_p);
        temp2=target_beta(id,TAU,lambda0,T,Z,X,Beta,N,C,gamma,k,Beta_curr,n_sample,x1_prev,uniqueid,no_p);
        lgratio=temp1-temp2;
        
        if (lgratio>=log(u))
        {
            for(int j=0; j<x1_curr.length(); j++){
                Beta_curr(k,j)=x1_curr[j];
            }
        }
        else
        {
            for(int j=0; j<x1_prev.length(); j++){
                Beta_curr(k,j)=x1_prev[j];
            }
                
        }
        return (Beta_curr(k,_));
}


// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]


double target_beta2(NumericVector id,NumericVector id2,double TAU,NumericVector lambda0,NumericVector T,
                    NumericVector T2,NumericVector Z,NumericVector Z2,NumericMatrix X,NumericMatrix X2,
                    NumericMatrix Beta,int N,int N2,NumericVector C,NumericVector C2,
                   double gamma,int k,NumericMatrix omega,NumericVector x1_new,double alpha,double n_sample,
                   NumericVector uniqueid,NumericVector uniqueid2,int no_p,int no_p2)
{
        double sum = 0.0, func=0.0;
        int m;
        NumericMatrix tempp;
        NumericMatrix temp=omega;
        for(m=1; m<n_sample; m++)
        {
                if(m==k) {
                        temp(m,_)=x1_new;
                }
                else
                {
                        temp(m,_)=omega(m,_);
                }
        }
        //int tempCols = std::extent<decltype(temp), 1>::value;
        for (int s=0; s<temp(m,_).size(); s++){
            tempp(s,0)=temp(m,s);
        }
        func=func+log(R::dgamma(temp(m,0),1,1,0));
        sum=LL2(id,TAU,lambda0,T,Z,X,tempp,N,C,gamma,alpha,uniqueid,no_p)+
                LL(id2,TAU,lambda0,T2,Z2,X2,tempp,N2,C2,gamma,uniqueid2,no_p2)+func;
        return (sum);
}


// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

NumericVector metroplis_beta2(NumericVector id,NumericVector id2,double TAU,NumericVector lambda0,NumericVector T,
                              NumericVector T2,NumericVector Z,NumericVector Z2,NumericMatrix X,
                              NumericMatrix X2,NumericMatrix Beta,int N,int N2,NumericVector C,NumericVector C2,
                             double gamma,int k,NumericMatrix Beta_curr,NumericVector x1_prev,double alpha,
                             double n_sample,double scale,NumericVector uniqueid,NumericVector uniqueid2,int no_p,int no_p2)
{
        double lgratio=0, u ,temp1=0, temp2=0;
        NumericVector x1_curr;
        u=R::runif(0,1);
        for(int j=0; j<x1_prev.length(); j++){
                x1_curr[j]=x1_prev[j]+scale*R::rnorm(0,1);  
        }
        
        temp1=target_beta2(id,id2,TAU,lambda0,T,T2,Z,Z2,X,X2,Beta,N,N2,C,C2,
                           gamma,k,Beta_curr,x1_curr,alpha,n_sample,uniqueid,uniqueid2,no_p,no_p2);
        temp2=target_beta2(id,id2,TAU,lambda0,T,T2,Z,Z2,X,X2,Beta,N,N2,C,C2,
                           gamma,k,Beta_curr,x1_prev,alpha,n_sample,uniqueid,uniqueid2,no_p,no_p2);
        lgratio=temp1-temp2;
        
        if (lgratio>=log(u))
        {
                Beta_curr(k,_)=x1_curr;
        }
        else
        {
                
                Beta_curr(k,_)=x1_prev;
                
        }
        return (Beta_curr(k,_));
}


// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

double target_lambda0(NumericVector id,NumericVector id2,double TAU,NumericMatrix lambda0,NumericVector T,
                      NumericVector T2,NumericVector Z,NumericVector Z2,NumericMatrix X,
                      NumericMatrix X2,NumericMatrix Beta,int N,int N2,NumericVector C,NumericVector C2,
                      double gamma,int k,NumericVector x1_new,double alpha,
                      double n_sample,NumericVector uniqueid,NumericVector uniqueid2,int no_p,int no_p2)
{
        double sum = 0.0, func=0.0;
        int m;
        NumericMatrix temp=lambda0;
        NumericVector tempp;
        for(m=1; m<n_sample; m++)
        {
                if(m==k) {
                        temp(m,_)=x1_new;
                }
                else
                {
                        temp(m,_)=lambda0(m,_);
                }
                
        }
        int tempCols = std::extent<decltype(temp), 1>::value;
        func=func+log(R::dgamma(temp(m,tempCols-1),1000,1000,0)); 
        for (int s=0; s<temp(m,_).size(); s++){
            tempp[s]=temp(m,s);
        }
        
        sum=LL2(id,TAU,tempp,T,Z,X,Beta,N,C,gamma,alpha,uniqueid,no_p)+
             LL(id2,TAU,tempp,T2,Z2,X2,Beta,N2,C2,gamma,uniqueid2,no_p2)+
             func;
        return(sum);
}


// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

NumericVector metroplis_lambda0(NumericVector id,NumericVector id2,double TAU,NumericMatrix lambda0,NumericVector T,
                         NumericVector T2,NumericVector Z,NumericVector Z2,NumericMatrix X,
                         NumericMatrix X2,NumericMatrix Beta,int N,int N2,NumericVector C,NumericVector C2,
                         double gamma,int k,NumericMatrix lambda0_curr,NumericVector x1_prev,double alpha,
                         double n_sample,double scale,NumericVector uniqueid,NumericVector uniqueid2,
                         int no_p,int no_p2)
        
{
        
        
        double lgratio=0, u ,temp1=0, temp2=0;
        NumericVector x1_curr;
        u=R::runif(0,1);
        for(int j=0; j<x1_prev.length(); j++){
                x1_curr[j]=x1_prev[j]+scale*R::rnorm(0,1);  
        }
        temp1=target_lambda0(id,id2,TAU,lambda0,T,T2,Z,Z2,X,X2,Beta,N,N2,C,C2,gamma,k,x1_curr,
                             alpha,n_sample,uniqueid,uniqueid2,no_p,no_p2);
        temp2=target_lambda0(id,id2,TAU,lambda0,T,T2,Z,Z2,X,X2,Beta,N,N2,C,C2,gamma,k,x1_prev,
                             alpha,n_sample,uniqueid,uniqueid2,no_p,no_p2);
        lgratio=temp1-temp2;
        
        if (lgratio>=log(u))
        {
                lambda0_curr(k,_)=x1_curr;
        }
        else
        {
                lambda0_curr(k,_)=x1_prev;
        }
        return(lambda0_curr(k,_));
}


// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]


double target_gamma(NumericVector id,double TAU,NumericVector lambda0,NumericVector T,
                    NumericVector Z,NumericMatrix X,NumericMatrix Beta,int N,NumericVector C,
                    NumericVector g,int k,double x1_new,double n_sample,NumericVector uniqueid,int no_p)
{
        double sum = 0.0, func=0.0;
        int m;
        NumericVector temp=g;
        for(m=1; m<n_sample; m++)
        {
                if(m==k) {
                        temp[m]=x1_new;
                }
                else
                {
                        temp[m]=g[m];
                }
        }
        
        func=func+log(R::dgamma(temp[m],1,1,0));
        sum=LL(id,TAU,lambda0,T,Z,X,Beta,N,C,
               temp[m],uniqueid,no_p)+func;
        return (sum);
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

double metroplis_gamma(NumericVector id,double TAU,NumericVector lambda0,NumericVector T,
                       NumericVector Z,NumericMatrix X,NumericMatrix Beta,int N,NumericVector C,
                       NumericVector g,int k,NumericVector gamma_curr,double x1_prev,
                       double n_sample,double scale,NumericVector uniqueid,int no_p)
{
        double lgratio,temp1,temp2;
        double x1_curr;
        double u=R::runif(0,1);
        double sd = 1.0;
        double mean = 0.0;
        double rnorm = normalRandom()*sd+mean;
        x1_curr=x1_prev+scale*rnorm;
        temp1=target_gamma(id,TAU,lambda0,T,Z,X,Beta,N,C,
                           g,k,x1_curr,n_sample,uniqueid,no_p);
        temp2=target_gamma(id,TAU,lambda0,T,Z,X,Beta,N,C,
                           g,k,x1_prev,n_sample,uniqueid,no_p);
        lgratio=temp1-temp2;
        
        if (lgratio>=log(u))
        {
                gamma_curr[k]=x1_curr;
        }
        else
        {
                gamma_curr[k]=x1_prev;
        }
        
        return  gamma_curr[k];
}


// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]


double target_gamma2(NumericVector id,NumericVector id2,double TAU,NumericVector lambda0,NumericVector T,
                     NumericVector T2,NumericVector Z,NumericVector Z2,NumericMatrix X,NumericMatrix X2,
                     NumericMatrix Beta,int N,int N2,NumericVector C,NumericVector C2,NumericVector g,int k,
                    double x1_new,double alpha,double n_sample,NumericVector uniqueid,NumericVector uniqueid2,
                    int no_p,int no_p2)
{
        double sum = 0.0, func=0.0;
        int m;
        NumericVector temp=g;
        for(m=1; m<n_sample; m++)
        {
                if(m==k) {
                        temp[m]=x1_new;
                }
                else
                {
                        temp[m]=g[m];
                }
        }
        
        func=func+log(R::dgamma(temp[m],1000,1000,0));
        sum=LL2(id,TAU,lambda0,T,Z,X,Beta,N,C,temp[m],alpha,uniqueid,no_p)+
             LL(id2,TAU,lambda0,T2,Z2,X2,Beta,N2,C2,temp[m],uniqueid2,no_p2)+func;
        return (sum);
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

double metroplis_gamma2(NumericVector id,NumericVector id2,double TAU,NumericVector lambda0,NumericVector T,
                        NumericVector T2,NumericVector Z,NumericVector Z2,NumericMatrix X,NumericMatrix X2,
                        NumericMatrix Beta,int N,int N2,NumericVector C,NumericVector C2,
                        NumericVector g,int k,NumericVector gamma_curr,double x1_prev,double alpha,
                        double n_sample,double scale,NumericVector uniqueid,NumericVector uniqueid2,
                        int no_p,int no_p2)
{
        double lgratio=0,temp1=0,temp2=0;
        double x1_curr;
        double u=R::runif(0,1);
        double sd = 1.0;
        double mean = 0.0;
        double rnorm = normalRandom()*sd+mean;
        x1_curr=x1_prev+scale*rnorm;
        temp1=target_gamma2(id,id2,TAU,lambda0,T,T2,Z,Z2,X,X2,
                            Beta,N,N2,C,C2,g,k,x1_curr,alpha,n_sample,uniqueid,uniqueid2,
                            no_p,no_p2);
        temp2=target_gamma2(id,id2,TAU,lambda0,T,T2,Z,Z2,X,X2,
                            Beta,N,N2,C,C2,g,k,x1_prev,alpha,n_sample,uniqueid,uniqueid2,
                            no_p,no_p2);
        lgratio=temp1-temp2;
        
        if (lgratio>=log(u))
        {
                gamma_curr[k]=x1_curr;
        }
        else
        {
                gamma_curr[k]=x1_prev;
        }
        
        return  gamma_curr[k];
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

List MH(NumericVector id,NumericVector TAU,NumericMatrix lambda0,NumericVector T,
        NumericVector Z,NumericMatrix X,NumericMatrix Beta,NumericMatrix Beta_sample,int N,NumericVector C,
        NumericVector gamma,double n_sample,NumericVector scale,NumericVector uniqueid,int no_p){
        int i;
        double u;
        double lgratio,temp1,temp2,temp3,x1_prev1,tau_pat;
        NumericVector x1_prev,temp,omega,omega_curr,gamma_curr,y_new;
        NumericMatrix Beta_curr,lambda0_curr,BETA;
        GetRNGstate();
        printf("%f\n",n_sample);
        for (i = 1; i < n_sample; i++)
        {
            Beta_curr=Beta_sample;
            for (int j =0;j<Beta_sample(i-1,_).size();j++){
                BETA(j,0)=Beta_sample(i-1,j);
            }
                //double sd = 1.0;
                //double mean = 0.0;
                //double rnorm = normalRandom()*sd+mean;
                //y_new=b1[i-1]+scale[0]*rnorm;
                
                printf("%d\n",i);
                u=R::runif(0,1);
                for (int j =0;j<lambda0(i-1,_).size();j++){
                    y_new[j]=lambda0(i-1,j)+scale[2]*R::rnorm(0,0.01);
                }
            temp1=LL(id,TAU[i-1],y_new,T,Z,X,BETA,N,C,gamma[i-1],uniqueid,no_p);
            printf("%f\n",temp1);
            temp2=LL(id,TAU[i-1],lambda0(i-1,_),T,Z,X,BETA,N,C,gamma[i-1],uniqueid,no_p);
            lgratio=temp1-temp2;
            if (lgratio>=log(u))  
            {for (int j =0;j<lambda0(i-1,_).size();j++){
                lambda0(i,j)=y_new[j];}
            }
            else {lambda0(i,_)=lambda0(i-1,_);}
            
               //Beta_curr=Beta_sample;
               printf("%d\n",i-1);
        for (int j =0;j<Beta_sample(i,_).size();j++){
               x1_prev[j]=Beta_sample(i-1,j);
               Beta_curr(i-1,j)=Beta_sample(i-1,j)+scale[4]*R::rnorm(0,1);
               }
        //    temp=metroplis_beta(id,TAU[i-1],lambda0(i,_),T,Z,X,BETA,N,C,gamma[i-1],i,
          //  Beta_curr,x1_prev,n_sample,scale[4],uniqueid,no_p);
              
            // for (int j =0;j<Beta_sample(i,_).size();j++){
              //   Beta_sample(i,j)=temp[j];
            //      y_new[j]=Beta(i-1,j)+scale[4]*R::rnorm(0,1);
              //}
              //Beta_curr.row(i)=y_new;
              //temp1=LL(id,TAU[i-1],lambda0(i,_),T,Z,X,Beta_curr,N,C,gamma[i-1],uniqueid,no_p);
              //printf("%f\n",temp1);
              //temp2=LL(id,TAU[i-1],lambda0(i,_),T,Z,X,Beta,N,C,gamma[i-1],uniqueid,no_p);
              //lgratio=temp1-temp2;
              //if (lgratio>=log(u))  
              //{Beta(i,_)=y_new;}
              //else {Beta(i,_)=Beta(i-1,_);}
                printf("%f\n",Beta_sample(i,0));
                printf("%f\n",Beta_sample(i,1));
            for (int j =0;j<Beta_sample(i,_).size();j++){
                    BETA(j,0)=Beta_sample(i,j);
                }
                
             omega_curr=NumericVector(N,R::rgamma(1,1));
                omega=omega_curr;
                for (int j = 0; j <no_p; j++)
                {
                    x1_prev1=omega[j];
                  printf("%f\n",x1_prev1);
                  //omega_curr[j]=x1_prev1+scale[4]*R::rnorm(0,0.01);
                    temp3=metroplis_omega(id,TAU[i-1],lambda0(i,_),T,Z,X,BETA,N,C,gamma[i-1],j,
                                          omega,x1_prev1,scale[4],uniqueid,no_p);
                    omega[j]=temp3;
                }
                //omega=center_dep(omega,no_p);
               
                gamma_curr=gamma;
                printf("%f\n",gamma_curr[i-1]);
                x1_prev1=gamma[i-1];
                printf("%f\n",x1_prev1);
                temp3=metroplis_gamma(id,TAU[i-1],lambda0(i,_),T,Z,X,BETA,N,C,gamma,i,
                gamma_curr,x1_prev1,n_sample,scale[4],uniqueid,no_p);
                printf("%f\n",temp3);
                gamma[i]=temp3;
                printf("%f\n", gamma[i]);
                
              
        tau_pat=target_tau(omega,no_p);
        TAU[i]=1/R::rgamma(0.001+N/2,1/(0.001+0.5*tau_pat));
        //TAU[i]=1/R::rgamma(0.001+N/2,1/(0.001+0.5));
        printf("%f\n", TAU[i]);
        }     
        
        PutRNGstate();
        
        return List::create(
                _["omega"]=omega,
                _["Beta"]= Beta_sample,
                _["tau"]= TAU,
                _["lambda0"]= lambda0,
                _["gamma"]= gamma);
        
}


// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

List MH2(NumericVector id,NumericVector id2,NumericVector TAU,NumericMatrix lambda0,NumericVector T,
         NumericVector T2,NumericVector Z,NumericVector Z2,NumericMatrix X,NumericMatrix X2,
         NumericMatrix Beta,NumericMatrix Beta_sample,int N,int N2,NumericVector C,NumericVector C2,
        NumericVector gamma,double n_sample,NumericVector alpha,NumericVector scale,NumericVector uniqueid,NumericVector uniqueid2,
        int no_p,int no_p2){
        int j;
        double u;
        double lgratio,temp1,temp2,temp3,x1_prev1,tau_pat;
        NumericVector x1_prev,temp,omega,omega_curr,gamma_curr,y_new,alpha_curr;
        NumericMatrix lambda0_curr,Beta_curr,BETA;;
        GetRNGstate();
        printf("%f\n",n_sample);
        for (int i = 1; i < n_sample; i++)
        {
                
                Beta_curr=Beta_sample;
            for (int j =0;j<Beta_sample(i-1,_).size();j++){
                BETA(j,0)=Beta_sample(i-1,j);
            }
                printf("%d\n",i);
                u=R::runif(0,1); 
                lambda0_curr=lambda0;
                for (int j =0;j<lambda0(i-1,_).size();j++){
                    x1_prev[j]=lambda0(i-1,j);
                    lambda0_curr(i-1,j)=lambda0(i-1,j)+scale[4]*R::rnorm(0,0.01);
                }
               
                //x1_prev=lambda0(i-1,_);
                temp=metroplis_lambda0(id,id2,TAU[i-1],lambda0,T,T2,Z,Z2,X,X2,BETA,
                 N,N2,C,C2,gamma[i-1],i,lambda0_curr,x1_prev,alpha[i-1],n_sample,scale[1],
                uniqueid,uniqueid2,no_p,no_p2);
                for (int j =0;j<lambda0(i,_).size();j++){  
                lambda0(i,j)=temp[j];
                } 
                
                
                
                //Beta_curr=Beta_sample;
                printf("%d\n",i-1);
                for (int j =0;j<Beta_sample(i,_).size();j++){
                    x1_prev[j]=Beta_sample(i-1,j);
                    Beta_curr(i-1,j)=Beta_sample(i-1,j)+scale[4]*R::rnorm(0,1);
                }
             //   temp=metroplis_beta2(id,id2,TAU[i-1],lambda0(i,_),T,T2,Z,Z2,X,X2,Beta,N,N2,C,C2,gamma[i-1],i,
            //Beta_curr,x1_prev,alpha[i-1],n_sample,scale[4],uniqueid,uniqueid2,no_p,no_p2);
              //  printf("%f\n",temp[0]);
                // for (int j =0;j<Beta_sample(i,_).size();j++){
                  //Beta_sample(i,j)=temp[j];
                  //  y_new[j]=Beta(i-1,j)+scale[4]*R::rnorm(0,1);
                //}
                //Beta_curr.row(i)=y_new;
                //temp1=LL(id,TAU[i-1],lambda0(i,_),T,Z,X,Beta_curr,N,C,gamma[i-1],uniqueid,no_p);
                //printf("%f\n",temp1);
                //temp2=LL(id,TAU[i-1],lambda0(i,_),T,Z,X,Beta,N,C,gamma[i-1],uniqueid,no_p);
                //lgratio=temp1-temp2;
                //if (lgratio>=log(u))  
                //{Beta(i,_)=y_new;}
                //else {Beta(i,_)=Beta(i-1,_);}
                printf("%f\n",Beta_sample(i,0));
                printf("%f\n",Beta_sample(i,1));
                for (int j =0;j<Beta_sample(i,_).size();j++){
                    BETA(j,0)=Beta_sample(i,j);
                }
                
                
              //  Beta_curr=Beta_sample;
                //printf("%d\n",i-1);
                //for (int j =0;j<Beta_sample(i-1,_).size();j++){
                  //  x1_prev[j]=Beta_sample(i-1,j);
                    //Beta_curr(i-1,j)=Beta_sample(i-1,j)+scale[4]*R::rnorm(0,1);
                //}
               
              //  temp=metroplis_beta2(id,id2,TAU[i-1],lambda0(i,_),T,T2,Z,Z2,X,X2,Beta,N,N2,C,C2,gamma[i-1],i,
            //                        Beta_curr,x1_prev,alpha[i-1],n_sample,scale[4],uniqueid,uniqueid2,no_p,no_p2);
              //  for (int j =0;j<Beta_sample(i,_).size();j++){
                //    Beta_sample(i,j)=temp[j];
                //}
               
               
                omega_curr=NumericVector(N,R::rgamma(1,1));
                omega=omega_curr;
                for (int j = 0; j <no_p; j++)
                {
                    x1_prev1=omega[j];
                    temp3=metroplis_omega(id,TAU[i-1],lambda0(i,_),T,Z,X,BETA,N,C,gamma[i-1],j,
                                          omega,x1_prev1,scale[4],uniqueid,no_p);
                   // temp3=metroplis_omega2(id,id2,TAU[i-1],lambda0(i,_),T,T2,Z,Z2,X,X2,BETA,N,N2,C,C2,gamma[i],j,omega_curr,x1_prev1,alpha[i-1],scale[4],uniqueid,uniqueid2,no_p,no_p2);
                    omega[j]=temp3;
                    printf("%f\n",omega[j]);
                }
                
                gamma_curr=gamma;
                
                x1_prev1=gamma[i-1];
                printf("%f\n",x1_prev1);
                //temp3=metroplis_gamma2(id,id2,TAU[i-1],lambda0(i,_),T,T2,Z,Z2,X,X2,BETA,N,N2,C,C2,
                 //                      gamma,i,gamma_curr,x1_prev1,alpha[i-1],n_sample,scale[4],uniqueid,uniqueid2,no_p,no_p2);
                 temp3=metroplis_gamma(id,TAU[i-1],lambda0(i,_),T,Z,X,BETA,N,C,gamma,i,
                                       gamma_curr,x1_prev1,n_sample,scale[4],uniqueid,no_p);
                
                gamma[i]=temp3;
                printf("%f\n",temp3);
               
               tau_pat=target_tau(omega,no_p);
                TAU[i]=1/R::rgamma(0.001+N/2,1/(0.001+0.5*tau_pat));
                printf("%f\n",TAU[i]);
                alpha_curr=NumericVector(N,R::runif(0,1));
                alpha=alpha_curr;
                for (int j = 0; j <no_p; j++){
                       x1_prev1=alpha[j];
                        
                temp3=metroplis_alpha(id,id2,TAU[i],lambda0(i,_),T,T2,Z,Z2,X,X2,BETA,N,N2,C,C2,
                       gamma[i],j,alpha_curr,x1_prev1,scale[4],uniqueid,uniqueid2,no_p,no_p2);
                        alpha[j]=temp3; 
                       printf("%f\n",alpha[j]);
                }
        }     
        
        PutRNGstate();
        
        return List::create(
                _["omega"]=omega,
                _["Beta"]= Beta_sample,
                _["tau"]= TAU,
                _["lambda0"]= lambda0,
                _["gamma"]= gamma,
                _["alpha"]= alpha);
        
}



// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

List MH3(NumericVector id,NumericVector id2,NumericVector TAU,NumericMatrix lambda0,NumericVector T,
         NumericVector T2,NumericVector Z,NumericVector Z2,NumericMatrix X,NumericMatrix X2,
         NumericMatrix Beta,NumericMatrix Beta_sample,int N,int N2,NumericVector C,NumericVector C2,
         NumericVector gamma,double n_sample,double alpha,NumericVector scale,
         NumericVector uniqueid,NumericVector uniqueid2,int no_p,int no_p2){
    int j;
    double u;
    double lgratio,temp1,temp2,temp3,x1_prev1,tau_pat;
    NumericVector x1_prev,temp,omega,omega_curr,gamma_curr,y_new,alpha_curr;
    NumericMatrix lambda0_curr,Beta_curr,BETA;;
    GetRNGstate();
    printf("%f\n",n_sample);
    for (int i = 1; i < n_sample; i++)
    {
        
        Beta_curr=Beta_sample;
        for (int j =0;j<Beta_sample(i-1,_).size();j++){
            BETA(j,0)=Beta_sample(i-1,j);
        }
        printf("%d\n",i);
        u=R::runif(0,1); 
        lambda0_curr=lambda0;
        for (int j =0;j<lambda0(i-1,_).size();j++){
            x1_prev[j]=lambda0(i-1,j);
            lambda0_curr(i-1,j)=lambda0(i-1,j)+scale[4]*R::rnorm(0,0.01);
        }
        
        //x1_prev=lambda0(i-1,_);
        temp=metroplis_lambda0(id,id2,TAU[i-1],lambda0,T,T2,Z,Z2,X,X2,BETA,
                               N,N2,C,C2,gamma[i-1],i,lambda0_curr,x1_prev,alpha,n_sample,scale[1],
                                                                                                    uniqueid,uniqueid2,no_p,no_p2);
        for (int j =0;j<lambda0(i,_).size();j++){  
            lambda0(i,j)=temp[j];
        } 
        
        
        
        //Beta_curr=Beta_sample;
        printf("%d\n",i-1);
        for (int j =0;j<Beta_sample(i,_).size();j++){
            x1_prev[j]=Beta_sample(i-1,j);
            Beta_curr(i-1,j)=Beta_sample(i-1,j)+scale[4]*R::rnorm(0,1);
        }
        //   temp=metroplis_beta2(id,id2,TAU[i-1],lambda0(i,_),T,T2,Z,Z2,X,X2,Beta,N,N2,C,C2,gamma[i-1],i,
        //Beta_curr,x1_prev,alpha[i-1],n_sample,scale[4],uniqueid,uniqueid2,no_p,no_p2);
        //  printf("%f\n",temp[0]);
        // for (int j =0;j<Beta_sample(i,_).size();j++){
        //Beta_sample(i,j)=temp[j];
        //  y_new[j]=Beta(i-1,j)+scale[4]*R::rnorm(0,1);
        //}
        //Beta_curr.row(i)=y_new;
        //temp1=LL(id,TAU[i-1],lambda0(i,_),T,Z,X,Beta_curr,N,C,gamma[i-1],uniqueid,no_p);
        //printf("%f\n",temp1);
        //temp2=LL(id,TAU[i-1],lambda0(i,_),T,Z,X,Beta,N,C,gamma[i-1],uniqueid,no_p);
        //lgratio=temp1-temp2;
        //if (lgratio>=log(u))  
        //{Beta(i,_)=y_new;}
        //else {Beta(i,_)=Beta(i-1,_);}
        printf("%f\n",Beta_sample(i,0));
        printf("%f\n",Beta_sample(i,1));
        for (int j =0;j<Beta_sample(i,_).size();j++){
            BETA(j,0)=Beta_sample(i,j);
        }
        
        
        //  Beta_curr=Beta_sample;
        //printf("%d\n",i-1);
        //for (int j =0;j<Beta_sample(i-1,_).size();j++){
        //  x1_prev[j]=Beta_sample(i-1,j);
        //Beta_curr(i-1,j)=Beta_sample(i-1,j)+scale[4]*R::rnorm(0,1);
        //}
        
        //  temp=metroplis_beta2(id,id2,TAU[i-1],lambda0(i,_),T,T2,Z,Z2,X,X2,Beta,N,N2,C,C2,gamma[i-1],i,
        //                        Beta_curr,x1_prev,alpha[i-1],n_sample,scale[4],uniqueid,uniqueid2,no_p,no_p2);
        //  for (int j =0;j<Beta_sample(i,_).size();j++){
        //    Beta_sample(i,j)=temp[j];
        //}
        
        
        omega_curr=NumericVector(N,R::rgamma(1,1));
        omega=omega_curr;
        for (int j = 0; j <no_p; j++)
        {
            x1_prev1=omega[j];
            temp3=metroplis_omega(id,TAU[i-1],lambda0(i,_),T,Z,X,BETA,N,C,gamma[i-1],j,
                                  omega,x1_prev1,scale[4],uniqueid,no_p);
            // temp3=metroplis_omega2(id,id2,TAU[i-1],lambda0(i,_),T,T2,Z,Z2,X,X2,BETA,N,N2,C,C2,gamma[i],j,omega_curr,x1_prev1,alpha[i-1],scale[4],uniqueid,uniqueid2,no_p,no_p2);
            omega[j]=temp3;
            printf("%f\n",omega[j]);
        }
        
        gamma_curr=gamma;
        
        x1_prev1=gamma[i-1];
        printf("%f\n",x1_prev1);
        //temp3=metroplis_gamma2(id,id2,TAU[i-1],lambda0(i,_),T,T2,Z,Z2,X,X2,BETA,N,N2,C,C2,
        //                      gamma,i,gamma_curr,x1_prev1,alpha[i-1],n_sample,scale[4],uniqueid,uniqueid2,no_p,no_p2);
        temp3=metroplis_gamma(id,TAU[i-1],lambda0(i,_),T,Z,X,BETA,N,C,gamma,i,
                              gamma_curr,x1_prev1,n_sample,scale[4],uniqueid,no_p);
        
        gamma[i]=temp3;
        printf("%f\n",temp3);
        
        tau_pat=target_tau(omega,no_p);
        TAU[i]=1/R::rgamma(0.001+N/2,1/(0.001+0.5*tau_pat));
        printf("%f\n",TAU[i]);
        
    }     
    
    PutRNGstate();
    
    return List::create(
        _["omega"]=omega,
        _["Beta"]= Beta_sample,
        _["tau"]= TAU,
        _["lambda0"]= lambda0,
        _["gamma"]= gamma);
    
}


