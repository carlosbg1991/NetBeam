clc; clear all; clear classes; close all;  %#ok
addpath('kriging/');  % Include Kriging folder (variogram and kriging)
set(0,'DefaultFigureColor','remove');  % No gray background in figures

%% SIMULATION CONFIGURATION
antIDList     = (1:1:4);  % Antenna ID, could be 1,2,3,4
expIDList     = (6:1:14);  % Experiment ID, could be 1,2,3,4,5
plotFlag      = false;  % Flag to plot results
N             = 12;  % Number of transmitter antennas
M             = 3;  % Number of receiver antennas
% Configuration for 1ST STAGE: DIRECT-VM 
config.n_samples     = 8;   % Maximum trials
config.new_additions = 1;   % New trials per iteration
config.newAddPoss    = 4;   % New possible possitions for DIRECT
config.nIter         = 1;  % Iterations to run Random over
% Configure comms environment
environment = 'outdoor';  % 'indoor', 'outdoor'
% Configure policy list
% policy_1stList = {'DIRECT-minVar','DIRECT-rand','UM','PI','random','DIRECT'};
policy_1stList = {'UM'};

%% Retreive data
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

H = zeros(M,N);  % Intermediate variable
H_exh = zeros(M,N);  % Output exhaustive

% Main loop
for idxPolicy = 1:length(policy_1stList)

    policy_1st = policy_1stList{idxPolicy};
    fprintf('Executing with policy "%s"\n',policy_1st);

    % Initialize Off-line learning stage sensibly...
    if strcmp(policy_1st,'DIRECT-minVar');      config.minSamples    = 2;   % Minimum Kriging initialization space
    elseif strcmp(policy_1st,'DIRECT-rand');    config.minSamples    = 4;   % Minimum Kriging initialization space
    elseif strcmp(policy_1st,'random');         config.minSamples    = 8;   % Minimum Kriging initialization space
    elseif strcmp(policy_1st,'PI');             config.minSamples    = 2;   % Minimum Kriging initialization space	
    elseif strcmp(policy_1st,'UM');             config.minSamples    = 4;   % Minimum Kriging initialization space
	else;  error('Policy not expected');
    end

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

    if strcmp(policy_1st,'random')
        save('CHANNEL_random','H');
    elseif strcmp(policy_1st,'UM')
        save('CHANNEL_UM','H');
    elseif strcmp(policy_1st,'EI')
        save('CHANNEL_EI','H');
    elseif strcmp(policy_1st,'PI')
        save('CHANNEL_PI','H');
    elseif strcmp(policy_1st,'DIRECT-rand')
        save('CHANNEL_DIRECT_rand','H');
    elseif strcmp(policy_1st,'DIRECT-minVar')
        save('CHANNEL_DIRECT_minVar','H');
    end
end

% Save exhaustive solution
H = H_exh;
save('CHANNEL_exhaustive','H');



% EOF