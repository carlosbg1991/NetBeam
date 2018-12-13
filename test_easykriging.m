%% Configure workspace
clear all; clear classes; close all; clc;  %#ok
addpath('kriging/');  % Include Kriging folder (variogram and kriging)
addpath('easyKriging/');  % Include Easy Kriging folder (simplified 1-dim)

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

K = 7;
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
grid minor
xlabel('Azimuth antenna steering','FontSize',12);
ylabel('Empirical channel gain','FontSize',12);
title('Kriging-based channel gain for FDBF','FontSize',12);
lg = legend([p1 p2 p3 p4],'Real measurement','Trials (known)','Prediction (mean)','Uncertainty (variance)');