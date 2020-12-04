function CBG_parse_channels(varargin)
% CBG_parse_channels - The script executes Kriging using the data from the
% experiments. It stores the channel matrices resulting from multiple
% acquisition functions (DIRECT-minVar, UM, PI, etc). Make sure to execute
% this script before assessing the global performance (test_global).
%
% Syntax:  CBG_parse_channels(environment)
%
% Inputs:
%    environment [optional] - String that takes the values 'indoor' and
%    'outdoor'. It selects the environment in which the experiments were
%    carried out. It defaults to indoor.
%
% Outputs: []
%
%
%------------- BEGIN CODE --------------



if (nargin==1)
    environment = varargin{1};
elseif (nargin==0)
    clear all; clear classes;  %#ok
    close all; clc;
    environment = 'outdoor';  % take indoor for example
else
    error('ERROR: The script only accepts 1 or none inputs\n');
end

fprintf('Selected environment: %s\n',environment);

%% Configure workspace
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
offlineList          = [2 4 6 8];  % Define the number of samples offline. 
                                   % Recall that DIRECT-UM consist of an
                                   % initial offline stage and then an 
                                   % advanced online stage. The more online                                   
                                   % samples, the more effective the 
                                   % algorithm is.
% Configure policy list
policy_1stList = {'DIRECT-minVar','DIRECT-rand','UM','PI','random'};
% policy_1stList = {'random'};

%% Retreive data
if ~exist('RESULTS','var')
    load('RESULTS','indoor','outdoor','paramList');
elseif ~exist('outdoor','var') || ~exist('indoor','var')
    CBG_parse_experiments;  % parse experimental DATA
end

%% Configure data
if strcmp(environment,'indoor')
    elevList    = indoor.elevList;
    azimList    = indoor.azimList;
    gainTot     = indoor.gainTot;
    chTot       = indoor.chTot;
elseif strcmp(environment,'outdoor')
    elevList    = outdoor.elevList;
    azimList    = outdoor.azimList;
    gainTot     = outdoor.gainTot;
    chTot       = outdoor.chTot;
end

% Main loop
for idxPolicy = 1:length(policy_1stList)

    policy_1st = policy_1stList{idxPolicy};
    fprintf('Policy: %s\n',policy_1st);
    
    for idxOffline = 1:length(offlineList)
        offline = offlineList(idxOffline);
        fprintf('\tOffline samples: %d\n',offline);
        
        % Configure offline trials
        config.minSamples = offline;
        
        H = zeros(M,N);  % Intermediate variable
        H_exh = zeros(M,N);  % Output exhaustive

        myAntIDList = [repmat((1:4),1,3) repmat((5:8),1,3) repmat((9:12),1,3)];
        myAntID = 1;
        for expID = expIDList
            for antID = antIDList
                Z = gainTot(:,:,antID,expID);  % Channel gain (exhaustive)

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
                H(RxID,myAntIDList(myAntID)) = chTot(idxAzim,idxElev,antID,expID);
                H_exh(RxID,myAntIDList(myAntID)) = chTot(idxAzim_exh,idxElev_exh,antID,expID);

                myAntID = myAntID + 1;  % Update indices
            end
        end

        save(sprintf('data/CHANNEL_%s_%s_%s-%s',environment,policy_1st,num2str(config.minSamples),num2str(config.n_samples)),'H');
    end
end

% Save exhaustive solution
H = H_exh;
save(sptrintf('CHANNEL_%s_exhaustive',environment),'H');



% EOF