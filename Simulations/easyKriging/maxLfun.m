function [thetaMLE,muMLE,sigmaMLE,lval] = maxLfun(samples,points,corFun,lowerTheta,upperTheta)
% maxlfun() returns the maximum likelihood values of mu, sigma and theta.
%
% Input:
% samples - nx by nSamples matrix
% points - spatial locations of samples
% corFun - desired correlation function: 'exp','sexp','poly','tri'
% lowerTheta (optional) - lower bound on correlation length
% upperTheta (optional) - upper bound on correlation length
%
% Output:
% thetaMLE - estimated correlation length
% muMLE - estimated mean
% sigmaMLE - estimated standard deviation
% lval - likelihood value
 
% Function to minimize using the Maximum Likelihood. theta is our
% variable. Covariance matrix given sample points, corr function and
% range
f = @(theta) -lfun(samples,calcCorrMat(points,corFun,theta));

% Compute the maximum likelihood by minimizing the distance between
% every point considering a particular correlation function. Returns
% the correlation range between samples
thetaMLE = fminbnd(f,lowerTheta,upperTheta);

% Compute the Likelihood of every point knowing their expected
% correlation range. Returns both, prediction and uncertainty.
[lval,muMLE,sigmaMLE] = lfun(samples,calcCorrMat(points,corFun,thetaMLE));

if thetaMLE == upperTheta
    printf('Consider increasing upperbound')
end
   
% EOF
