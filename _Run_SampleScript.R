rm(list=ls())

source("_LIB_Data_Simulator.R");
source("_Funcs_l1ML.R")
source("_LIB_Evaluator.R");

p1 = 30; p2 = 60; n = 100;

B = GenMat_Coef(p=p1,q=p2,sparsity=5/p1,mg_low=1,mg_high=2);
InvCov_X = GenMat_invcov(q=p1,sparsity = 5/p1,type = 'bidiag',target_condition_number = 5,mg_high = 0.5, mg_low = NULL);
InvCov_E = GenMat_invcov(q=p2,sparsity = 5/p2,type = 'bidiag',target_condition_number = 5,mg_high = 0.5, mg_low = NULL);

Data = GenData_2layer(n,InvCov_X=InvCov_X,B,InvCov_E,standardize_cov=TRUE,target_SNR=1.8);

lambda.opt = 0.123; rho.opt = 0.234;

# Only the first two arguments are necessary.
# Strongly recommend set lambda and rho by yourself.
result = l1ML_Main(Data$Y,Data$X,lambda=lambda.opt,rho=rho.opt,initializer="Lasso",screening=T,alpha=0.1,ss=T,nboot=20);
   
## Evaluation for B
Eval_sparse(true = Data$B, est = result$B.est, directed = TRUE)                     

## Evaluation for Theta
Eval_sparse(true = Data$InvCov_E, est = result$Theta.est, directed = FALSE)                     
