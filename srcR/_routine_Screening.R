# ******************************************
# Description: return the skeleton of the directed edges (in binary)
#				based on de-biased Lasso and simultaneous testing; 
# File type: routine to be sourced by _Funcs_l1ML.R
#
# Last modified on 02/01/2017
# *******************************************

source("_estimate_SSLasso.R");

ScreeningFunc = function(Ymat,Xmat,alpha=0.1,adjust=TRUE,method="BH")
{
    # Args:
    #(data) Ymat: response matrix, n by p2;
	#(data) Xmat: covariate matrix, n by p1;
	#(param) alpha: simultaneous testing cut-off level, default is 0.1;
	#(param) adjust: if simultaneous testing is performed, default is TRUE;
	#(param) method: correction method for multiple comparison, choose between "BH" or "bonferroni";
	# -----------------------
	# Returns:
    # a p1 by p2 binary matrix indicating if the edge should be included;
	
	p1 = ncol(Xmat); p2 = ncol(Ymat); n = nrow(Xmat);
	
	if (nrow(Ymat)!=n)
    {
		stop("Error in Screening(): sample sizes for responses and covariates differ.\n");
	}
	
	p.val = array(0,c(p1,p2))
	output.list_prior = foreach(j = 1:p2, .export = c('SSLasso', 'adjust', 'Lasso', 'slim', 'NoiseSd', 'InverseLinfty', 'InverseLinftyOneRow', 'SoftThreshold')) %dopar%
    {
			f = SSLasso(Xmat,Ymat[,j],verbose=FALSE,intercept=FALSE);
			p.val_j = f$pvals;
			p.val_j # the p-values of the j-th regression
    }
    
	for (j in 1:p2)
    {
		p.val[,j] = output.list_prior[[j]]; 
	}
	
	if (adjust)
		p.val = array(p.adjust(as.vector(p.val),method=method),c(p1,p2));
	
	skeleton = (p.val <= alpha);
	return(skeleton);
}
