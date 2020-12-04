% TEST_ORIENTATION - The script loads results from a wide variety of
% experiments (indoors and outdoors) containing the channel measurements
% for the all the combinations of azimuth and elevation angles. Then, it
% enforces baseline antenna orientation policies (RANDOM, PI, UM) and
% compares them to the proposed antenna orientation policy in NetBeam
% (DIRECT-UM or DIRECT-minVar). The results show the gap to optimality in
% terms of channel gain.
%
% This script generates Fig. 11 of the publication :
% [1] C. Bocanegra, K. Alemdar, S. Garcia, C. Singhal and K. R. Chowdhury,
%     “NetBeam: Network of Distributed Full-dimension
%     Beamforming SDRs for Multi-user Heterogeneous Traffic,” IEEE Dynamic
%     Spectrum (DySpan), Newark, NJ, 2019
%
% Syntax:  test_orientation([])
%
% Inputs: []
%
% Outputs: []
%
%
%------------- BEGIN CODE --------------


%% Configure workspace
clc; clear all; clear classes; close all;  %#ok
addpath('../../Beamforming/');  % SDP solver in Beamforming repository
addpath('kriging/');  % SDP solver in Beamforming repository
addpath('data/');  % where results are stored (to be loaded) 
set(0,'DefaultFigureColor','remove');  % No gray background in figures

%% SIMULATION CONFIGURATION
antIDList            = (1:4);  % Antenna ID, could be 1,2,3,4
expIDList            = (6:10);  % Experiment ID, could be 1,2,3,4,5
plotFlag             = false;  % Flag to plot results
N                    = 12;  % Number of transmitter antennas
M                    = 3;  % Number of receiver antennas
SNRdemands           = [0.7; 0.6; 0.8];  % Minimum SINR for each user
% Configuration for 1ST STAGE: DIRECT-VM 
config.n_samples     = 8;   % Maximum trials
config.new_additions = 2;   % New trials per iteration
config.newAddPoss    = 4;   % New possible possitions for DIRECT
config.nIter         = 10;  % Iterations to run Random over
offlineList          = [2 4 6 8];  % Define the number of samples offline. 
                                   % Recall that DIRECT-UM consist of an
                                   % initial offline stage and then an 
                                   % advanced online stage. The more online                                   
                                   % samples, the more effective the 
                                   % algorithm is.
% Configure comms environment
environment          = 'indoor';  % 'indoor', 'outdoor'
% Define antenna orientation policies
policy_1stList       = {'random','DIRECT-rand','DIRECT-minVar','PI','UM'};

%% PARSE Data if not done before
if ~exist('RESULTS','var')
    load('RESULTS','indoor','outdoor','paramList');
elseif ~exist('outdoor','var') || ~exist('indoor','var')
    CBG_parse_experiments;  % parse experimental DATA
end

%% Configure data
if strcmp(environment,'indoor')
    DATA = indoor;
elseif strcmp(environment,'outdoor')
    DATA = outdoor;
end
elevList    = DATA.elevList;  % Local copy
azimList    = DATA.azimList;  % Local copy

gap2OptAv = zeros(length(offlineList),length(policy_1stList));
gap2OptMax = zeros(length(offlineList),length(policy_1stList));
gap2OptSum = zeros(length(offlineList),length(policy_1stList));
for idxPolicy1 = 1:length(policy_1stList)
    policy_1st = policy_1stList{idxPolicy1};
    fprintf('Policy: %s\n',policy_1st);
    for idxOffline = 1:length(offlineList)
        offline = offlineList(idxOffline);
        % Configure offline trials
        config.minSamples = offline;
        %% 1 stage: Antenna orientation
        myAntIDList = [repmat((1:4),1,3) repmat((5:8),1,3) repmat((9:12),1,3)];
        myAntID = 1;
        for expID = expIDList
            for antID = antIDList
                Z = DATA.gainTot(:,:,antID,expID);  % Channel gain (exhaustive)

                %% 1 stage: Maximize antenna gain per tx-rx pair
                [pred_gain,pred_elevation,pred_azimuth,~,~,~,...
                 Zhat,Zvar,Zexh] = CBG_kriging(Z,elevList,azimList,config,policy_1st,plotFlag);

                % find optimum angles - Exhaustive
                [exhv_gain, max_idx] = max(Z(:));
                [idx_exhv_maxGain_X,idx_exhv_maxGain_Y] = ind2sub(size(Z),max_idx);
                [X,Y] = meshgrid(elevList,azimList); 
                exhv_elevation = X(idx_exhv_maxGain_X,idx_exhv_maxGain_Y);
                exhv_azimuth = Y(idx_exhv_maxGain_X,idx_exhv_maxGain_Y);

                %% 1.5 stage: Retrieve channel information from selected angular configuration
                [~,idxElev] = min(abs(pred_elevation-elevList));
                [~,idxAzim] = min(abs(pred_azimuth-azimList));
                [~,idxElev_exh] = min(abs(exhv_elevation-elevList));
                [~,idxAzim_exh] = min(abs(exhv_azimuth-azimList));

                % Store final results in channel matrix
                myExpID = expID - 5;  % first 5 are dummy data
                RxID = mod(myExpID-1,3) + 1;
                H(RxID,myAntIDList(myAntID)) = DATA.chTot(idxAzim,idxElev,antID,expID);
                H_exh(RxID,myAntIDList(myAntID)) = DATA.chTot(idxAzim_exh,idxElev_exh,antID,expID);

                myAntID = myAntID + 1;  % Update indices
            end
        end

        % Determine gap to optimum
        gap2OptAv(idxOffline,idxPolicy1) = mean(abs(H_exh) - abs(H),'all');
        gap2OptMax(idxOffline,idxPolicy1) = max(abs(H_exh) - abs(H),[],'all');
        gap2OptSum(idxOffline,idxPolicy1) = sum(abs(H_exh) - abs(H),'all');
    end
end
%% AVERAGE
groupLabels = {'Random', 'DIRECT-RD', 'DIRECT-UM', 'PI', 'UM'};
stackData = gap2OptAv.';
plotBarStackGroups(stackData, groupLabels,1);
lg = legend('Off-line: 2 trials / On-line: 8 trials',...
            'Off-line: 4 trials / On-line: 8 trials',...
            'Off-line: 6 trials / On-line: 8 trials',...
            'Off-line: 8 trials');
set(lg,'FontSize',8);
ylabel('(Mean) Gap to optimality (linear gain)','FontSize',12);
grid minor;
% Modify colors
a = findobj(gca,'type','bar');
a(4).FaceColor = [0 191 255]./255;  %light green
a(3).FaceColor = [70 130 180]./255;  %dark green
a(2).FaceColor = [106 90 205]./255;  %light blue
a(1).FaceColor = [0 0 128]./255;  %dark blue
%% MAX
groupLabels = {'Random', 'DIRECT-RD', 'DIRECT-UM', 'PI', 'UM'};
stackData = gap2OptMax.';
plotBarStackGroups(stackData, groupLabels,2);
lg = legend('Off-line: 2 trials / On-line: 8 trials',...
            'Off-line: 4 trials / On-line: 8 trials',...
            'Off-line: 6 trials / On-line: 8 trials',...
            'Off-line: 8 trials');
set(lg,'FontSize',8);
ylabel('(Max) Gap to optimality (linear gain)','FontSize',12);
grid minor;
%% SUM
groupLabels = {'Random', 'DIRECT-RD', 'DIRECT-UM', 'PI', 'UM'};
stackData = gap2OptSum.';
plotBarStackGroups(stackData, groupLabels,3);
lg = legend('Off-line: 2 trials / On-line: 8 trials',...
            'Off-line: 4 trials / On-line: 8 trials',...
            'Off-line: 6 trials / On-line: 8 trials',...
            'Off-line: 8 trials');
set(lg,'FontSize',8);
ylabel('(Sum) Gap to optimality (linear gain)','FontSize',12);
grid minor;