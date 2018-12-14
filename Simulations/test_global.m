%% Configure workspace
clc; clear all; clear classes; close all;  %#ok
addpath('../../Beamforming/');  % SDP solver in Beamforming repository
set(0,'DefaultFigureColor','remove');  % No gray background in figures

%% SIMULATION CONFIGURATION
antIDList            = (1:1:4);  % Antenna ID, could be 1,2,3,4
expIDList            = (6:1:14);  % Experiment ID, could be 1,2,3,4,5
plotFlag             = true;  % Flag to plot results
N                    = 12;  % Number of transmitter antennas
M                    = 3;  % Number of receiver antennas
SNRdemands           = [0.7; 0.6; 0.8];  % Minimum SINR for each user
% Configuration for 1ST STAGE: DIRECT-VM 
offline              = 2;  % Number of offline trials
% Configuration for 2ND STAGE: Antenna selection
% Configuration for 3RD STAGE: SDB beamforming
sigma2               = ones(12,1);  % Noise variance
Pt_max               = 1;  % Maximum transmitted power per each radius
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


policy_1stList = {'random','DIRECT-rand_2','DIRECT-rand_4','DIRECT-minVar_8','PI_4','PI_2','DIRECT-minVar_4','DIRECT-minVar_2','UM_4','UM_2'};
legends        = {'Random','DIRECT-rand_2','DIRECT-rand_4','DIRECT_8','PI_4','PI_2','DIRECT-minVar_4','DIRECT-minVar_2','UM_4','UM_2'};
policy_2ndList = {'optimum','greedy-equtv','random'};

txPowerTot = [];
rxSNRsum = [];
for idxPolicy2 = 1:length(policy_2ndList)
    txPowerTot1 = [];
    rxSNRsum1 = [];
    for idxPolicy1 = 1:length(policy_1stList)
        policy_1st = policy_1stList{idxPolicy1};
        policy_2nd = policy_2ndList{idxPolicy2};

        %% 1 stage: Antenna orientation (pre-stored channel)
        % To get the new channels, please run "CBG_parse_channels.m"
        if strcmp(policy_1st,'DIRECT-minVar_2')
            ws = load('CHANNEL_DIRECT_minVar_2-8.mat');
            H = ws.H;
        elseif strcmp(policy_1st,'DIRECT-minVar_4')
            ws = load('CHANNEL_DIRECT_minVar_2-8.mat');
            H = ws.H;
        elseif strcmp(policy_1st,'DIRECT-minVar_8')
            ws = load('CHANNEL_DIRECT_minVar_8.mat');
            H = ws.H;
        elseif strcmp(policy_1st,'DIRECT-rand_2')
            ws = load('CHANNEL_DIRECT_rand_2-8.mat');
            H = ws.H;
        elseif strcmp(policy_1st,'DIRECT-rand_4')
            ws = load('CHANNEL_DIRECT_rand_4-8.mat');
            H = ws.H;
        elseif strcmp(policy_1st,'random')
            ws = load('CHANNEL_random_8.mat');
            H = ws.H;
        elseif strcmp(policy_1st,'PI_2')
            ws = load('CHANNEL_PI_2-8.mat');
            H = ws.H;
        elseif strcmp(policy_1st,'PI_4')
            ws = load('CHANNEL_PI_4-8.mat');
            H = ws.H;
        elseif strcmp(policy_1st,'UM_2')
            ws = load('CHANNEL_UM_2-8.mat');
            H = ws.H;
        elseif strcmp(policy_1st,'UM_4')
            ws = load('CHANNEL_UM_4-8.mat');
            H = ws.H;
        elseif strcmp(policy_1st,'exhaustive')
            ws = load('CHANNEL_exhaustive.mat');
            H = ws.H;
        else
            error('Policy not expected');
        end

        %% 2 stage: Antenna selection
        [finalAssign,~,assignation] = CBG_antSel(H,SNRdemands,policy_2nd);

        %% 3 stage: Beamforming
        repeat = true;
        SNRdemands1 = SNRdemands;
        while repeat
            [w, SNR, repeat] = CBG_sdp_solver(H, N, M, Pt_max, SNRdemands1, sigma2, assignation);
            % Adjust demanded SNR and try it again
            [maxSNR,i] = max(SNRdemands1);
            minSNR = min(SNRdemands1);
            delta = 0.05;
            SNRdemands1(i) = SNRdemands1(i) - delta;
        end
        
        %% Final: show results
        fprintf('Policy 1: %s, Policy 2: %s -> SNR: | ',policy_1st,policy_2nd);
        fprintf('%.7f |',SNR);
        fprintf('\n');
        fprintf('Power used: %.5f\n',sum(sum(abs(w).^2)));
        
        %% Append results - Tx power
        txPowerTot1 = [txPowerTot1 sum(sum(abs(w).^2))];  %#ok<AGROW>
        
        [~,lostID] = find(SNR - SNRdemands.' > 1e-4);
        if isempty(lostID); lostSNR = 0;
        else;               lostSNR = SNRdemands(lostID).' - SNR(lostID);
        end
        rxSNRsum1 = [rxSNRsum1 sum(lostSNR)];  %#ok<AGROW>
    end
    txPowerTot = [txPowerTot ; txPowerTot1];  %#ok<AGROW>
    rxSNRsum = [rxSNRsum ; rxSNRsum1];  %#ok<AGROW>
end

%% PLOT OVERALL TRANSMIT POWER
groupLabels = policy_2ndList;
stackData = txPowerTot;
plotBarStackGroups(stackData, groupLabels,1);
lg = legend(legends);
set(lg,'FontSize',8);
title('analysis on the MAX GAIN achieved','FontSize',12);
ylabel('Overall transmit power (linear)','FontSize',12);
grid minor;
% Modify colors
a = findobj(gca,'type','bar');
a(1).FaceColor = [178 255 102]./255;  %light green
a(2).FaceColor = [0 153 76]./255;  %dark green
a(3).FaceColor = [51 153 255]./255;  %light blue
a(4).FaceColor = [0 0 204]./255;  %dark blue
a(5).FaceColor = [255 255 102]./255;  %light yellow
a(6).FaceColor = [153 153 0]./255;  %dark yellow
a(7).FaceColor = [255 102 102]./255;  %light red
a(8).FaceColor = [255 0 0]./255;  %red
a(9).FaceColor = [153 0 0]./255;  %dark red
a(10).FaceColor = [0 0 0]./255;  %black

groupLabels = policy_2ndList;
stackData = rxSNRsum;
plotBarStackGroups(stackData, groupLabels,2);
lg = legend(legends);
set(lg,'FontSize',8);
title('analysis on the MAX GAIN achieved','FontSize',12);
ylabel('Overall SNR lost','FontSize',12);
grid minor;
% Modify colors
a = findobj(gca,'type','bar');
a(1).FaceColor = [178 255 102]./255;  %light green
a(2).FaceColor = [0 153 76]./255;  %dark green
a(3).FaceColor = [51 153 255]./255;  %light blue
a(4).FaceColor = [0 0 204]./255;  %dark blue
a(5).FaceColor = [255 255 102]./255;  %light yellow
a(6).FaceColor = [153 153 0]./255;  %dark yellow
a(7).FaceColor = [255 102 102]./255;  %light red
a(8).FaceColor = [255 0 0]./255;  %red
a(9).FaceColor = [153 0 0]./255;  %dark red
a(10).FaceColor = [0 0 0]./255;  %black