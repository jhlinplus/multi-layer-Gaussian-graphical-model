
rm(list=ls())

source("utilsR/data_simulator.R");
source("utilsR/evaluator.R");
source("srcR/_estimate_l1ML.R")

###############################
#### Generate Some Synthetic Data
###############################

p1 = 30; p2 = 60; n = 100;

B = gen_sparse_coef_mtx(p=p1,q=p2,sparsity=5/p1,mg_low=1,mg_high=2);
InvCov_X = gen_sparse_inv_cov(q=p1,sparsity = 5/p1,type = 'bidiag',target_condition_number = 5,mg_high = 0.5, mg_low = NULL);
InvCov_E = gen_sparse_inv_cov(q=p2,sparsity = 5/p2,type = 'bidiag',target_condition_number = 5,mg_high = 0.5, mg_low = NULL);

myData = generate_data_2layerGGM(n,InvCov_X=InvCov_X,B,InvCov_E,standardize_cov=TRUE,target_SNR=1.8);

###############################
#### Perform Estimation
###############################

lambda.opt = 0.123; rho.opt = 0.234;

# Only the first two arguments are necessary.
# Strongly recommend set lambda and rho by yourself.
result = l1ML_main(myData$Y,myData$X,lambda=lambda.opt,rho=rho.opt,initializer="Lasso",screening=T,alpha=0.1,ss=T,nboot=20);

###############################
#### Perform Evaluation
###############################

## Evaluation for B
eval_sparse(true = myData$B, est = result$B.est, directed = TRUE)

## Evaluation for Theta
eval_sparse(true = myData$InvCov_E, est = result$Theta.est, directed = FALSE)
