function [krige,CIupper,CIlower] = krigeIt(condPoints,condVals,uncondPoints,corFun,mu,sigma,theta,CIalpha)
% krigeIt returns kriging interpolation and confidence intervals using
% multivariate Gaussian maximum likelihood estimation to estimate the mean 
% and autocorrelation from the data
%
% Mandatory inputs:
% condPoints (nx by 1 vector) - spatial or time location of data
% condVals (nx by 1 vector) - value of process at those points
% uncondPoints (nx by 1 vector) - spatial or time location where prediction
%                                 is desired
% corFun - correlation function. Options: 'exp','sexp','poly','tri'
% mu - mean of process underlying data (if known)
% sigma - standard deviation of process underlying data (if known)
% theta - correlation length of process underlying data (if known)

% Optional inputs:
% CIalpha - confidence level for confidence intervals (default: 0.05)
%
% Outputs:
% krige - kriging interpolation of data
% CIupper - upper confidence interval
% CIlower - lower confidence interval
% errorVar - kriging variance

if nargin < 8
    CIalpha = 0.05;
end
    
pointsExt = [uncondPoints; condPoints];

nCond = length(condPoints);
nUncond = length(uncondPoints);
condInd = [nUncond+1:nUncond+nCond];

corrMat = calcCorrMat(pointsExt,corFun,theta);
covMat = sigma^2*corrMat;
C = covMat(condInd,condInd);
b = covMat(condInd,:);
beta = C\b;
krigeExt = mu + beta'*(condVals - mu);

krige = krigeExt(1:nUncond);

errorVarExt = sigma^2 - sum(b.*beta)';
errorVar = errorVarExt(1:nUncond);

CIupper = norminv(CIalpha/2,krige,sqrt(errorVar));
CIlower = norminv(1-CIalpha/2,krige,sqrt(errorVar));

end

