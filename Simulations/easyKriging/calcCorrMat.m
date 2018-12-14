function [ corMat ] = calcCorrMat( points,corFun, corLen )
%CALCCORRMAT Summary of this function goes here
%   Detailed explanation goes here

distVec = pdist(points);
adjDist = distVec/corLen;

switch corFun
    case 'exp'
        corVec = exp(-2*adjDist);
    case 'sexp'
        corVec = exp(-pi*adjDist.^2);
    case 'poly'
        corVec = (1 + adjDist).^-3;
    case 'tri'
        corVec = 1 - adjDist;
        corVec(corVec < 0) = 0;
end

corMat = squareform(corVec);
[m,~] = size(corMat);
corMat(1:m+1:end) = 1;

end

