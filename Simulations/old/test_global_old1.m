%% Configure workspace
clear all; clear classes; close all; clc;  %#ok
addpath('kriging/');  % Include Kriging folder (variogram and kriging)
set(0,'DefaultFigureColor','remove');  % No gray background in figures

%% PARSE Data if not done before
if ~exist('RESULTS','var')
    load('RESULTS','indoor','outdoor','paramList');
elseif ~exist('outdoor','var') || ~exist('indoor','var')
    CBG_parse_experiments;  % parse experimental DATA
end
elevList    = outdoor.elevList;  % Local copy
azimList    = outdoor.azimList;  % Local copy

%% SIMULATION CONFIGURATION
antIDList     = (1:1:4);  % Antenna ID, could be 1,2,3,4
expIDList     = (6:1:14);  % Experiment ID, could be 1,2,3,4,5
plotFlag      = false;  % Flag to plot results
N             = 12;  % Number of transmitter antennas
M             = 3;  % Number of receiver antennas
SNRdemands    = [0.1; 0.3; 0.15];  % Minimum SINR for each user
% Configuration for 1ST STAGE: DIRECT-VM 
config.n_samples     = 8;  % Maximum trials
config.minSamples    = 8;  % Minimum Kriging initialization space
config.new_additions = 1;  % New trials per iteration
config.newAddPoss    = 4;  % New possible possitions for DIRECT
config.nIter         = 1;  % Iterations to run Random over
policy_1st = 'DIRECT-minVar';
% Configuration for 2ND STAGE: Antenna selection
policy_2nd    = 'optimum';
% Configuration for 3RD STAGE: SDB beamforming
sigma2        = ones(12,1);  % Noise variance
Pt_max        = 1;  % Maximum transmitted power per each radius

H = zeros(M,N);
myAntIDList = [repmat((1:4),1,3) repmat((5:8),1,3) repmat((9:12),1,3)];
myAntID = 1;
for expID = expIDList
    for antID = antIDList
        Z = outdoor.gainTot(:,:,antID,expID);  % Channel gain (exhaustive)

        %% 1 stage: Maximize antenna gain per tx-rx pair
        [pred_gain,pred_elevation,pred_azimuth,...
         exhv_gain,exhv_elevation,exhv_azimuth,...
         Zhat,Zvar,Zexh] = CBG_kriging(Z,elevList,azimList,config,policy_1st,plotFlag);

        %% 1.5 stage: Retrieve channel information from selected angular configuration
        [~,idxElev] = min(abs(exhv_elevation-elevList));
        [~,idxAzim] = min(abs(exhv_azimuth-azimList));
        
        % Store final results in channel matrix
        myExpID = expID - 5;
        RxID = mod(myExpID-1,3) + 1;
        H(RxID,myAntIDList(myAntID)) = outdoor.chTot(idxAzim,idxElev,antID,expID);
        
        myAntID = myAntID + 1;  % Update indices
    end
end

%% 2 stage: Antenna selection
[finalAssign,~,assignation] = CBG_antSel(H,SNRdemands,policy_2nd);

%% 3 stage: Beamforming
[w, SNR] = CBG_sdp_solver(H, N, M, Pt_max, SNRdemands, sigma2, assignation);