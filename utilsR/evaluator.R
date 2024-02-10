Eval_sparse = function(true,est,directed=FALSE)
{
    if (directed==FALSE)
    {
        true.v=upperTriangle(!!true, diag=FALSE);
        est.v=upperTriangle(!!est, diag=FALSE);
    }
    else
    {
        true.v=as.vector(!!true); 
        est.v=as.vector(!!est);
    }
    
    TP = sum(est.v & true.v); 
    TN = sum((!est.v) & (!true.v));
    FP = sum(est.v & !true.v); 
    FN = sum(!est.v & true.v);
    
    SEN = TP/(TP+FN); # aka true positive rates
    SPC = TN/(TN+FP); # aka true negative rate;
    DE_MCC = sqrt((TP+FP))*sqrt((TP+FN))*sqrt((TN+FP))*sqrt((TN+FN));
    MCC = (TP*TN-FP*FN)/DE_MCC;
    Ferror = norm(true-est,"F")/norm(true,"F");
    
    result = c(SEN,SPC,MCC,Ferror);
    names(result) = c("SEN","SPC","MCC","ERROR")
    
    return(result)
}