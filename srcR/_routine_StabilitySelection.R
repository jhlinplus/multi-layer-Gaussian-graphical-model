# ***********************************************************************
# Description: perform stability selection on precision matrix, 
#				only probability matrix is return, no cut-off is set.
#
# Reference: "Stability Selection", N.Meinshausen; P.Buhlmann, JRSSB 2010; 
#
# File type: routine to be sourced by _Funcs_l1ML.R;
#
# Last modified on 02/01/2017 by JL.
# ***********************************************************************

pkgs = c("huge","gdata");
new = pkgs[!(pkgs%in%installed.packages()[,"Package"])]
if (length(new)){
    for (pkg in new)
    {
        install.packages(pkg, dependencies = TRUE)
    }
}
sapply(pkgs, require, character.only = TRUE);

StabilitySelection = function(X,method="glasso",rho.seq=NULL,nbootstrap=50,VERBOSE=TRUE)
{
	# Args:
	# (data) X: n*p centered data matrix for inverse covariance matrix estimation;
	# (param) method: underlying method for estimating the skeleton, choose between mb and glasso;
	# (param) rho.seq: the sequence of tuning parameters; default is (0.2 to 3)*sqrt(log(p)/n);
	# (param) nbootstrap: number of bootstrap samples; 
	# (param) VERBOSE: if the program is noisy ...
	# ----------------------
	# Returns: 
	# a p*p symmetric probability matrix, each entry is the existing probability of an edge.
	
	if (VERBOSE)
	{
		cat("Stability selection with a total number of", nbootstrap, "bootstrapped samples ...")
	}
	
	n = nrow(X); p = ncol(X);
	halfsize = as.integer(n/2); # each bootstrap sample has half the sample size.
	
	if (is.null(rho.seq))
		rho.seq = seq(from=0.2,to=3,length=20)*sqrt(log(p)/n);
	
	freq = matrix(0,length(rho.seq),p*(p-1)/2);
	if (method=="mb")
	{
		for (i in 1:(nbootstrap/2))
		{
			perm = sample(n);
			i1 = perm[1:halfsize];
			i2 = perm[(halfsize+1):n];
		
			Theta.est = huge(X[i1,],lambda=rho.seq,method=method,verbose=FALSE)$path;
			for (j in 1:length(rho.seq))
			{
				freq[j,] = freq[j,] + upperTriangle(as.matrix(Theta.est[[j]]));
			}
		
			Theta.est = huge(X[i2,],lambda=rho.seq,method=method,verbose=FALSE)$path;
			for (j in 1:length(rho.seq))
			{
				freq[j,] = freq[j,] + upperTriangle(as.matrix(Theta.est[[j]]));
			}
		
			if (VERBOSE){
				if (i==round(nbootstrap/8)) cat("25%...");
				if (i==round(nbootstrap/4)) cat("50%...");
				if (i==round(3*nbootstrap/8)) cat("75%...");
				if (i==nbootstrap/2) cat("Done!\n");
			}
		}
	}
	else if (method=="glasso")
    {
		for (i in 1:(nbootstrap/2))
        {
            perm = sample(n);
			i1 = perm[1:halfsize];
			i2 = perm[(halfsize+1):n];
			
            ## estimate on the first half of the sample
			glasso.est = huge(X[i1,],lambda=rho.seq,method="glasso",verbose=FALSE)$icov;
			
            for (j in 1:length(rho.seq))
            {
				freq[j,] = freq[j,] + upperTriangle(abs(sign(as.matrix(glasso.est[[j]]))));
			}
		
            ## estimate on the second half of the sample
			glasso.est = huge(X[i2,],lambda=rho.seq,method="glasso",verbose=FALSE)$icov;
			
            for (j in 1:length(rho.seq))
            {
				freq[j,] = freq[j,] + upperTriangle(abs(sign(as.matrix(glasso.est[[j]]))));
			}
            
            if (VERBOSE){
                if (i == round(nbootstrap/8)) cat("25%...");
                if (i == round(nbootstrap/4)) cat("50%...");
                if (i == round(3*nbootstrap/8)) cat("75%...");
                if (i == nbootstrap/2) cat("Done!\n");
            }
		}
	}	
		
	freq = freq/nbootstrap;
	
    result = array(0,c(p,p))
    upperTriangle(result) = apply(freq,2,max)
    result = result + t(result)
    diag(result) = 1
    
    return(result)
}
