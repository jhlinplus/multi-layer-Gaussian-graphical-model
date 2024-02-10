# ************************************
# Description: collections of routines that are necessary for l1ML_main();
#
# File type: Script for routines; 
# Routines: 
#	Calculate.Obj() -- calculate the objective function;
#	Calculate.BIC() -- calculate the BIC value;
#
# Last modified on 02/01/2017;
# ************************************

# require(gdata);

Calculate.Obj = function(Y,X,B,Theta,lambda,rho)
{
	n = nrow(X); p1 = ncol(X); p2 = ncol(Y);
	
	loss = -log(det(Theta)) + sum( diag(t(Y-X%*%B) %*% (Y-X%*%B) %*% Theta) )/n;
	penalize.B = lambda*sum(abs(B));
	penalize.Theta = rho*(sum(abs(Theta))-sum(diag(Theta)));
	
	objval = loss + penalize.B + penalize.Theta;
	return(objval);
}

Calculate.BIC = function(X,Y,B,Theta)
{
	# use KL loss; 
	# BIC = -log(det(Theta)) + tr(S%*%Theta) + sum(diag(t(Y-X%*%B) %*% (Y-X%*%B) %*% Theta))/n + (log(n)/n)*(nonzeros in upperTriangle(Theta) + nonzeros in B)
	# we prefer smaller BIC value
	
	n = nrow(X);
	KL.loss = sum(diag(t(Y-X%*%B) %*% (Y-X%*%B) %*% Theta))/n - log(det(Theta)) + sum(diag(t(Y-X%*%B) %*% (Y-X%*%B) %*% Theta))/n;
	penalization = log(n)/n*(sum(abs(upperTriangle(Theta))>1e-6) + sum(abs(B)>1e-6));
	
	bicval = KL.loss + penalization;
	return(bicval); 
}
