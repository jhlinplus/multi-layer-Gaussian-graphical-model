# ******************************************
# Description: library of functions that generates the data in
# "Penalized Maximum Likelihood Estimation of Multi-layered Gaussian Graphical Models"
# J.Lin; S.Basu; M.Banerjee; G.Michailidis; 2016 JMLR
#
# Last modified on 07/28/2019 by JL
# *******************************************

pkgs = c("matrixcalc","MASS","gdata");
new = pkgs[!(pkgs%in%installed.packages()[,"Package"])]
if (length(new)){
    for (pkg in new)
    {
        install.packages(pkg, dependencies = TRUE)
    }
}
sapply(pkgs, require, character.only = TRUE);

# generates a p*q regression coefficient matrix for model y = B'x + e
# where y is q times 1, B is p times q, x is p times 1
GenMat_Coef = function(p,q,sparsity=NULL,mg_low,mg_high){
    
	# p regressors, q regressions
    # default sparsity is 5/p
	if (is.null(sparsity))
		sparsity = 5/p;
	
	B = array(0,c(p,q));
    for (i in 1:p)
        for (j in 1:q)
            B[i,j] = rbinom(1,1,sparsity)*sample(c(-1,1),1)*runif(1,mg_low,mg_high);
	
	return(B)
}

# GenMat_invcov() generates an inverse covariance matrix
GenMat_invcov = function(q,sparsity = NULL,type = 'random',target_condition_number = 5,mg_high = 1, mg_low = 0.5)
{
    Omega = array(0,c(q,q));
    diag(Omega) = 0;
    
    
    if (type == "random")
    {
        if (is.null(sparsity))
            sparsity = 2/q;
        
        if (is.null(mg_low) || is.null(mg_high))
            stop("mg_high or mg_low cannot be NULL when type is set to random.\n");
        
        for (i in 1:(q-1))
        {
            for (j in (i+1):q)
            {
                Omega[j,i] = rbinom(1,1,sparsity)*sample(c(-1,1),1)*runif(1,mg_low,mg_high);
                Omega[i,j] = Omega[j,i];
            }
        }
    }
    else if (type == "bidiag")
    {
        if (is.null(mg_high))
            stop("mg_high cannot be NULL when type is set to bidiag.\n");
        
        for (i in 1:(q-1)){
            for (j in (i+1):q)
            {
                Omega[j,i] = ifelse(j-i==1,mg_high,0);
                Omega[i,j] = Omega[j,i];
            }
        }
    }
    else if (type == "tridiag")
    {
        if (is.null(mg_low) || is.null(mg_high))
            stop("mg_high or mg_low cannot be NULL when type is set to tridiag.\n");
        
        for (i in 1:(q-1)){
            for (j in (i+1):q){
                Omega[j,i] = ifelse(j-i==1,mg_high,ifelse(j-i==2,mg_low,0));
                Omega[i,j] = Omega[j,i];
            }
        }
    }
    
    Omega = Omega + diag(q)*(0.001+abs(min(eigen(Omega)$values)));
    egval = eigen(Omega)$values;
    CN_Omega = max(egval)/min(egval)
    
    while( CN_Omega > target_condition_number )
    {
        Omega = Omega + 0.005*diag(q);
        egval = eigen(Omega)$values;
        CN_Omega = max(egval)/min(egval);
    }
    
    return(Omega)
}

GenData_2layer = function(n,InvCov_X=NULL,B,InvCov_E,standardize_cov=TRUE,target_SNR=NULL){
    
    ## generate a model of the following form:
    ## Y = XB + E;
    # Args:
    # (@param) n: sample size
    # (@param) InvCov_X: inverse covariance matrix of X; default to diag(p)
    # (@param) B: X -> Y regression coefficient matrix
    # (@param) InvCov_E:inverse covariance matrix of E; cannot be NULL
    # (@param) standardize: if we standardize the covariance matrix, default=TRUE
    # (@param) target_SNR: if specified, will rescale B to match the target SNR. Set with caution.
    # >>>>>
    # Returns:
    # list(), layer 1 data X, layer 2 data Y; model parameters after possible scaling.
    
    p1 = nrow(B); p2 = nrow(InvCov_E);
    
    SigmaE = solve(InvCov_E);
    
    if (is.null(InvCov_X))
        SigmaX = diag(p1)
    else
        SigmaX = solve(InvCov_X);
    
    if (standardize_cov)
    {
		SigmaE = cov2cor(SigmaE);
        InvCov_E = solve(SigmaE);
        InvCov_E = ifelse(abs(InvCov_E)>1e-6,InvCov_E,0);
        
        SigmaX = cov2cor(SigmaX);
        InvCov_X = solve(SigmaX);
        InvCov_X = ifelse(abs(InvCov_X)>1e-6,InvCov_X,0);
    }
	
	X = mvrnorm(n,mu=rep(0,p1),Sigma=SigmaX);
	E = mvrnorm(n,mu=rep(0,p2),Sigma=SigmaE);
	
	# generate Y
    signal = X %*% B;
	Y = signal + E;
	
    if ( !is.null(target_SNR) )
    {
        current_SNR = mean( apply(signal,2,sd)/apply(E,2,sd) );
        cat(sprintf("SNR before scaling = %.2f.\n", current_SNR));
        
        B = B * target_SNR/current_SNR;
        signal = X %*% B;
        Y = signal + E;
    }
    cat(sprintf("Final SNR  = %.2f.\n", mean(apply(signal,2,sd)/apply(E,2,sd)) ));
    
    return( list(X=X,Y=Y,B=B,InvCov_X=InvCov_X,InvCov_E=InvCov_E) );
}




























