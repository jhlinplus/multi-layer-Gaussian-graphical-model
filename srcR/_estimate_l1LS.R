# ********************************************************
# Description: estimate directed edges using penalized least square and 
#				undirected edges using glasso based on residuals;
#
# File type: Script for main function;
# Required files:
#	_Routine_Screening.R;
#   _Funcs_SSLasso.R;
#
# References: http://www.jmlr.org/papers/volume17/16-004/16-004.pdf
# Last modified on 02/01/2019 by JL.
#
# ********************************************************

pkgs = c('glmnet','huge');
new = pkgs[!(pkgs%in%installed.packages()[,"Package"])]
if (length(new)){
    for (pkg in new)
    {
        install.packages(pkg, dependencies = TRUE)
    }
}
sapply(pkgs, require, character.only = TRUE);

source("_routine_Screening.R");

l1LS_main = function(Y,X,lambda=NULL,rho=NULL,initializer="Lasso",skeleton=NULL)
{
	# Argsv:
	# (data) Y: response matrix, n by p2;
	# (data) X: covariate matrix, n by p1;
	# (param) lambda: tuning parameter for the regularization term, default 0.2*sqrt(log(p1)/n);
	# (param) rho: tuning parameter for glasso; default sqrt(log(p2)/n);
	# (param) initializer: penalized LS method, choose between "Lasso" and "Ridge";
    # (param) skeleton: skeleton of the regression coefficients; if NULL then no skeleton is provided.
	# -----------------------
	# Returns (list):
	# B0: estimated regression coefficients matrix;
	# Theta0: estimated undirected edges in the second layer; 
	
	n = nrow(X); p1 = ncol(X); p2 = ncol(Y); 
	
	if (is.null(lambda))
    {
		lambda = 0.2*sqrt(log(p1)/n);
	}
	if (is.null(rho))
    {
		rho = sqrt(log(p2)/n);
	}
	
	if (!is.null(skeleton))
    {
        Existing.edges = diag(1:p1) %*% skeleton;
    }
	else
	{
		cat("No skeleton is detected.\n");
		skeleton = array(1,c(p1,p2));
		Existing.edges = diag(1:p1) %*% skeleton;
	}
	
	# estimate through Lasso or Ridge
	if (initializer=="Lasso")
	{
		output.list_LS = foreach(j=1:p2, .export = c('glmnet')) %dopar% 
		{
			B_j = rep(0,p1);
			if (length(which(Existing.edges[,j]!=0))==0)
			{
				B_j = rep(0,p1);
			}	
			else if (length(which(Existing.edges[,j]!=0))==1)
			{
				B_j[which(Existing.edges[,j]!=0)] = lm(Y[,j]~X[,which(Existing.edges[,j]!=0)]+0)$coef;
			}	
			else{	
				temp = glmnet(X[,which(Existing.edges[,j]!=0)],Y[,j],intercept=FALSE);
				B_j[which(Existing.edges[,j]!=0)] = predict(temp,s=lambda,type="coefficients")[-1];
			}
			B_j
		}
	}
	else if (initializer=="Ridge")
	{
		output.list_LS = foreach(j=1:p2)%dopar%
		{
			B_j = rep(0,p1)
			if (length(which(Existing.edges[,j]!=0))==0)
			{
				B_j = rep(0,p1);
			}
			else if (length(which(Existing.edges[,j]!=0))==1)
			{
				B_j[which(Existing.edges[,j]!=0)] = lm(Y[,j]~X[,which(Existing.edges[,j]!=0)]+0)$coef;
			}
			else
			{	
				temp = glmnet(X[,which(Existing.edges[,j]!=0)],Y[,j],intercept=FALSE,alpha=0);
				B_j[which(Existing.edges[,j]!=0)] = predict(temp,s=lambda,type="coefficients")[-1];
			}
			B_j
		}
	}
	else
	{
		stop("Error in l1LS_Main(): invalid method, choose between Lasso and Ridge.\n");
	}
	
	## collect the result
	B.est = array(0,c(p1,p2));
	for (j in 1:p2)
	{
		B.est[,j] = output.list_LS[[j]]
	}
	
    ## estimate inverse covariance matrix using glasso based on residuals
	ResMat = Y - X %*% B.est;
	Theta.est = as.matrix(huge(ResMat,rho,method="glasso",verbose=FALSE)$icov[[1]]);
	
	return(list(B0 = B.est,Theta0 = Theta.est))
}
