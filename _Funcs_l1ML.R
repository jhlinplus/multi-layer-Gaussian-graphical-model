# ********************************************************
# Description: estimate undirected edges within layers and directed edges between layers
# "Penalized Maximum Likelihood Estimation of Multi-layered Gaussian Graphical Models"  
# J.Lin; S.Basu; M.Banerjee; G.Michailidis; 2016 JMLR 
# 
# File type: main function;
#
# Required files (in the working directory):
#
#   _Funcs_l1LS.R;
#	_Routine_screening.R;
#   _Funcs_SSLasso.R;
#	_Routine_StabilitySelection.R;
#   _Routine_Aux.R;
#
# Last modified on 02/01/2017 by JL.
# ********************************************************

pkgs = c("glmnet","glasso","doParallel");
new = pkgs[!(pkgs%in%installed.packages()[,"Package"])]
if (length(new)){
    for (pkg in new)
    {
        install.packages(pkg, dependencies = TRUE)
    }
}
sapply(pkgs, require, character.only = TRUE);
registerDoParallel(cores=detectCores(all.tests = FALSE, logical = TRUE)-1);
cat(sprintf("Total number of workers = %d.\n",getDoParWorkers()));

source("_Funcs_l1LS.R");
source("_Routine_StabilitySelection.R");
source("_Routine_Aux.R");

l1ML_Main = function(Y,X,lambda=NULL,rho=NULL,initializer="Lasso",screening=T,alpha=0.1,method.correction="BH",ss=T,method.ss="glasso",rho.ss=NULL,nboot=50,VERBOSE=FALSE,iter_max=100,manualbreak=50)
{
    # Args:
    # (data) Y: data in layer 2, n by p2;
    # (data) X: data in layer 1, n by p1;
    # (param) lambda: lambda_n in eqn(2), default sqrt(log(p1)/n), use at your own discretion!
    # (param) rho: rho_n in eqn(2), default sqrt(log(p2)/n), use at your own discretion!
    # (param) initializer: type of penalized LS to use for initialization, choose between "Lasso" and "Ridge";
    # (param) screening: indicator if the screening procedure is included, default is TRUE;
    # (param) alpha: cut-off for p-values, when screening is on.
    # (param) method.correction: method used for multiple comparison, choose between "BH" and "bonferroni";
    # (param) ss: indicator if the stability selection is performed, default is TRUE;
    # (param) method.ss: underlying method to use for stability selection, choose between "glasso" and "mb";
    # (param) rho.ss: the sequence of tuning parameters on which stability selection is performed;
    # (param) nboot: number of bootstrapped samples;
    # (param) VERBOSE: if we track the change of the objective function, default = FALSE;
    # (param) iter_max: maximum iteration allowed for the update between B and Theta;
    # (param) manualbreak: specify iteration no. to break the alternate update; if NULL, wait until convergence;
    # --------------------------------------
    # Returns (list):
    # B.est: estimated directed edges between layer 1 and 2, p1 by p2 matrix.
    # Theta.est: estimated undirected edges within layer 2, p2 by p2 matrix.
    # BICvalue: BIC value with the current results plugged-in.

    n = nrow(X); p1 = ncol(X); p2 = ncol(Y);
    
    if (nrow(Y)!=n)
    {
        stop("Error in l1ML_Main(): sample sizes in Layer 1 and Layer 2 diff.\n")
    }
    
    if (is.null(lambda))
    {
        lambda = sqrt(log(p1)/n);
    }
    if (is.null(rho))
    {
        rho = sqrt(log(p2)/n);
    }
    ######################
    ### Screening
    ######################
    if (screening){
        cat("Step 0: screening is on, proceed with debiased Lasso in conjunction with",method.correction,"correction, the cut-off for p-values is set at",as.character(alpha),".\n");
        skeleton = ScreeningFunc(Ymat=Y,Xmat=X,alpha=alpha,adjust=TRUE,method=method.correction);
        
    }
    else{
        cat("Step 0: Screening is off.\n");
        skeleton = array(1,c(p1,p2));
    }
    Existing.edges = diag(1:p1) %*% skeleton; # intermediate quantity for subsequent use;
    
    ######################
    ### Initialization
    ######################
    cat("Step 1: Penalized LS initialization with",initializer,".\n");
    LS = l1LS_Main(Y,X,lambda=NULL,rho=NULL,initializer=initializer,skeleton=skeleton);
    B.initial = LS$B0; Theta.initial = LS$Theta0;
    
    ###############################
    ### Prior to alternate update
    ###############################
    Obj.initial = Calculate.Obj(Y,X,B.initial,Theta.initial,lambda,rho);
    Objval = c(); iter = 0; CONVERGE=F;
    updateTheta=FALSE; # suppress Theta update until B is relatively stable.
    
    # residual matrix
    Rmat = Y - X %*% B.initial;
    # r_j's on p.10, expression in between eqn(7) and (8);
    R = array(0,c(n,p2));
    for (j in 1:p2){
        for (k in (1:p2)[-j]){
            R[,j] = R[,j] + Theta.initial[j,k]*Rmat[,k];
        }
        R[,j] = R[,j]/Theta.initial[j,j];
    }
    ###############################
    ### Alternate update
    ###############################
    cat("Step 2: Alternate update.\n");
    B_new = B.initial; Theta_new = Theta.initial;
    
    while(!CONVERGE)
    {
        iter = iter + 1;
        if (iter > iter_max){
            stop("Error in l1ML_Main(): iteration overflow for during the alternate update.\n");
        }
        if (!is.null(manualbreak)){
            if (iter > manualbreak){
                warning("l1ML_Main(): alternate update is broken manually at iteration =",manualbreak,"\n");
                break;
            }
        }
        
        B_old = B_new; Theta_old = Theta_new;
        
        # update B
        output.list_B = foreach (j=1:p2)%dopar%
        {
            B_j = rep(0,p1);
            if (length(which(Existing.edges[,j]!=0))==0){
                B_j = B_j;
            }
            else if (length(which(Existing.edges[,j]!=0))==1){
                coef1 = lm((Y[,j]+R[,j])~X[,which(Existing.edges[,j]!=0)]+0)$coef;
                ## this is the least square solution, need to apply a soft-thresholding
                coef1_st = max(0,coef1-lambda/(2*Theta_old[j,j]))*sign(coef1);
                B_j[which(Existing.edges[,j]!=0)] = coef1_st;
            }
            else if (length(which(Existing.edges[,j]!=0))>=2){
                temp = glmnet(X[,which(Existing.edges[,j]!=0)],(Y[,j]+R[,j]),intercept=FALSE);
                B_j[which(Existing.edges[,j]!=0)] = predict(temp,s=lambda/(2*Theta_old[j,j]),type="coefficients")[-1];
            }
            B_j;
        }
        for (j in 1:p2){
            B_new[,j] = output.list_B[[j]];
        }
        
        # Update Theta
        if (iter >=10 | norm(B_new-B_old,"F")<0.1 | updateTheta == TRUE){
            updateTheta = TRUE;
            ResMat = Y - X%*%B_new;
            Theta_new = as.matrix(huge(ResMat,rho,method="glasso",verbose=FALSE)$icov[[1]]);
        }
        
        Objval[iter] = Calculate.Obj(Y,X,B_new,Theta_new,lambda,rho);
        
        ## check convergence
        Obj_diff = ifelse(iter==1,Objval[1]-Obj.initial,Objval[iter] - Objval[iter-1])
        CONVERGE = (abs(Obj_diff)<1e-4);
        
        if (VERBOSE){
            cat(paste(c('iter=','Obj_diff='),c(iter,round(Obj_diff,5)),sep=" "),'\n');
        }
        
        
        Rmat = Y - X %*% B_new;
        R = array(0,c(n,p2));
        # update R_i's
        for (i in 1:p2){
            for (k in (1:p2)[-i]){
                R[,i] = R[,i]+ Theta_new[i,k]*Rmat[,k];
            }
            R[,i] = R[,i]/(Theta_new[i,i]);
        }
    }
    ############
    ### Refit
    ############
    cat("Step 3: Refitting B Matrix.\n");
    B.refit = array(0,c(p1,p2));
    output.list_refit = foreach (j = 1:p2) %dopar% {
        B.refit_j = rep(0,p1);
        if (length(which(abs(B_new[,j])>1e-6))>0){
            B.refit_j[which(abs(B_new[,j])>1e-6)]=lm((Y[,j]+R[,j])~X[,which(abs(B_new[,j])>1e-6)]+0)$coef;
        }
        B.refit_j;
    }
    for (j in 1:p2){
        B.refit[,j] = output.list_refit[[j]];
    }
    ResMat.refit = Y - X %*% B.refit;
    
    if (ss){
        cat("Step 4: ");
        if (is.null(rho.ss)){
            rho.ss = seq(from=0.2,to=3.0,length=10)*sqrt(log(p2)/n)
        }
        ProbMat = StabilitySelection(ResMat.refit,method=method.ss,rho.seq=rho.ss,nbootstrap=nboot,VERBOSE=TRUE);
        WeightMat = 1 - ProbMat;
        cat("Step 5: Refitting with weighted glasso ...");
        Theta.refit = glasso(var(ResMat.refit),2*rho*WeightMat,penalize.diagonal=FALSE)$wi;
        cat("Done.\n");
    }
    else if (!ss){
        cat("Step 4: Refitting with glasso.");
        Theta.refit = as.matrix(huge(ResMat.refit,rho,method="glasso",verbose=FALSE)$icov[[1]])
    }
    
    B.est = B.refit; Theta.est = Theta.refit;
    BICvalue = Calculate.BIC(X,Y,B.est,Theta.est);
    
    cat("Returns: B.est, Theta.est, BICvalue.\n");
    return(list(B.est = B.est,
                Theta.est = Theta.est,
                BICvalue = BICvalue,
                skeleton = skeleton,
                B.initial = B.initial,
                Theta.initial = Theta.initial,
                Objval = c(Obj.initial,Objval)
            ))
}
