%% Configure workspace
clc; clear all; clear classes; close all;  %#ok
addpath('../../Beamforming/');  % SDP solver in Beamforming repository
set(0,'DefaultFigureColor','remove');  % No gray background in figures

%% SIMULATION CONFIGURATION
antIDList            = (1:1:4);  % Antenna ID, could be 1,2,3,4
expIDList            = (6:1:14);  % Experiment ID, could be 1,2,3,4,5
plotFlag             = false;  % Flag to plot results
N                    = 12;  % Number of transmitter antennas
M                    = 3;  % Number of receiver antennas
SNRdemands           = [0.7; 0.6; 0.8];  % Minimum SINR for each user
% Configuration for 1ST STAGE: DIRECT-VM 
config.n_samples     = 16;   % Maximum trials
config.new_additions = 1;   % New trials per iteration
config.newAddPoss    = 4;   % New possible possitions for DIRECT
config.nIter         = 1;  % Iterations to run Random over
% Configure comms environment
environment          = 'outdoor';  % 'indoor', 'outdoor'

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


policy_1stList = {'random','DIRECT-rand','DIRECT-minVar','PI','UM'};

offlineList = [2 4 6 8];
count = 1;
for idxPolicy1 = 1:length(policy_1stList)
    policy_1st = policy_1stList{idxPolicy1};
    for offline = offlineList
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
        gap2Opt(count) = abs(H_exh) - abs(H);
        count = count + 1;  % counter for indexing results
    end
end