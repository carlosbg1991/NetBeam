%% Configure workspace
clear all; clear classes; close all; clc;  %#ok
addpath('kriging/');  % Include Kriging folder (variogram and kriging)

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

%% Call 1st stage - Kriging
num_trials = 30;
plotFlag = true;
[opt_krig,opt_exhv,Z_krigMean,Z_krigVari] = CBG_kriging(Z,elevList,azimList,num_trials,plotFlag);

%% Call 2nd stage - Antenna selection
K = nTxAntennas;  % Max number of Tx to be allocated to each Rx.
K = 2*floor(N/M);  % Let's assume we can allocate all the radios to one receiver
AntAlloc = CBG_antSel(chMax,K);

%% Call 3rd stage - Beamforming