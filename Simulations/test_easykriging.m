function test_easykriging(varargin)
% TEST_EASYKRIGING - The script loads results containing the channel
% measurements for the all the combinations of azimuth and elevation
% angles. Then, it showcases the benefits of krigging where a number of
% samples are taken as knowns and the model predicts the channel for
% unvisited angles, as well as the uncertainty associated to such
% prediction.
%
% This script generates Fig. 6 of the publication (make sure to call this
% function with knownRatio=7):
% [1] C. Bocanegra, K. Alemdar, S. Garcia, C. Singhal and K. R. Chowdhury,
%     “NetBeam: Network of Distributed Full-dimension
%     Beamforming SDRs for Multi-user Heterogeneous Traffic,” IEEE Dynamic
%     Spectrum (DySpan), Newark, NJ, 2019
%
% Syntax:  test_easykriging(knownRatio)
%
% Inputs:
%    knownRatio [optional] - Possitive integer that captures the ratio of
%    known angles over the unknowns. This should be positive and greater
%    than 1, i.e., 2 or greater. The lower knownRatio is, the lower the
%    uncertainty in the system is and the more accurate the predictions
%    are.
%
% Outputs: []
%
%
%------------- BEGIN CODE --------------

if (nargin==1)
    K = varargin{1};
elseif (nargin==0)
    clear all; clear classes; close all; clc;  %#ok
    K = 7;
else
    error('ERROR: The script only accepts 1 or none inputs\n');
end

%% Configure workspace
addpath('kriging/');  % Include Kriging folder (variogram and kriging)
addpath('easyKriging/');  % Include Easy Kriging folder (simplified 1-dim)
addpath('data/');  % Include Easy Kriging folder (simplified 1-dim)

%% Load experimental data
load('results-27-Nov-2018_01:32:43.mat','azimList','elevList','chTot','nTxAntennas','appliedElev','appliedAzym');
N = nTxAntennas;  % Number of transmitter antennas
M = 1;  % Number of receiver antennas

%% Parse it and create matrix with channel gain per azimuth-elevation pairs
chGain = abs(chTot);
[X,Y] = meshgrid(elevList,azimList);
dim = size(X,1)*size(X,2);
Z_prel = zeros(nTxAntennas,length(elevList)*length(azimList));
Z = zeros(length(azimList),length(elevList),nTxAntennas);
for id = 1:nTxAntennas
    for t = 1:dim
        elev = X(t);
        azym = Y(t);
        idx_elev = find(appliedElev==elev);
        idx_azym = find(appliedAzym==azym);
        idx = intersect(idx_elev,idx_azym);
        Z_prel(id,t) = mean(chGain(id,idx));
    end
    Z(:,:,id) = reshape(Z_prel(id,:),[size(X,1),size(X,2)]);
end

condPoints = Y(1:K:end,8);  % Sample Points
condVals = Z(1:K:end,8);  % Sample values

%% Apply Kriging and get prediction and confidence intervals for sampled data

lx = max(azimList);  % Max value to predict
nx = 1e3;  % Maximum number of points to predict
dx = lx/nx;  % Distance between prediction points
uncondPoints = [dx/2:dx:lx]';  % Interpolated range

corFun = 'sexp';  % Kernel type: 'exp','sexp','poly','tri'
lowerTheta = 0;  % Lower correlation bound (range for samples to be corr)
upperTheta = 100;  % Controls the overfitting (higher, smoother)

% Maximum Likelihood estimation
[theta,mu,sigma,lval] = maxLfun(condVals,condPoints,corFun,lowerTheta,upperTheta);

% Interpolate over the unknown values using Kriging. Returns predictions
% and confidence bounds
[krige,CIupper,CIlower] = krigeIt(condPoints,condVals,uncondPoints,corFun,mu,sigma,theta);

%% Plotting section
figure;
hold on
A = [uncondPoints; flipud(uncondPoints)];
B = [CIlower; flipud(CIupper)];
p4 = fill(A',B',[0.8,0.8,0.8]);
set(p4,'EdgeColor','none')
p1 = plot(Y(:,8),Z(:,8),'color',[0 0 1],'LineWidth',2);
p2 = scatter(condPoints,condVals,'MarkerEdgeColor',[1 0 0],'LineWidth',2,'Marker','o','SizeData',50,'LineWidth',2,'MarkerFaceColor',[1,1,1]);
p3 = plot(uncondPoints,krige,'color',[1 0 0],'LineWidth',2);
grid on;
xlabel('Azimuth antenna steering','FontSize',10);
ylabel('Empirical channel gain','FontSize',10);
title(sprintf('Kriging-based channel gain for FDBF, K=%d',K),'FontSize',10);
lg = legend([p1 p2 p3 p4],'Real measurement','Trials (known)','Prediction (mean)','Uncertainty (variance)');