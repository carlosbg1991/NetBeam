%% Configure workspace
clc; clear all; clear classes; close all;  %#ok
addpath('../../Beamforming/');  % SDP solver in Beamforming repository
set(0,'DefaultFigureColor','remove');  % No gray background in figures

%% SIMULATION CONFIGURATION
N                    = 12;  % Number of transmitter antennas
M                    = 3;  % Number of receiver antennas
SNRdemands           = [0.6; 0.4; 0.7];  % Minimum SINR for each user
% Configuration for 1ST STAGE: DIRECT-VM 
% Pre-stored
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

policy_1stList = {'DIRECT-minVar_4','UM_4','PI_4','DIRECT-minVar_8','DIRECT-rand_4','random'};
legends        = {'DIRECT-UM','UM','PI','DIRECT','DIRECT-RD','Random'};
policy_2ndList = {'optimum','greedy','random'};

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
title('INDOORS','FontSize',12);
ylabel('Overall tx. power (linear)','FontSize',12);
% Modify colors
a = findobj(gca,'type','bar');
a(6).FaceColor = [0 104 255]./255;
a(1).FaceColor = [76 76 76]./255;
a(2).FaceColor = [127 127 127]./255;
a(3).FaceColor = [178 178 178]./255;
a(4).FaceColor = [204 204 204]./255;
a(5).FaceColor = [229 229 229]./255;
xlim([0.5 3.5]);
pos = get(gcf, 'Position');
set(gcf,'position',[pos(1),pos(2),677,235]);
grid minor;

groupLabels = policy_2ndList;
stackData = rxSNRsum;
plotBarStackGroups(stackData, groupLabels,2);
lg = legend(legends);
set(lg,'FontSize',8);
title('INDOORS','FontSize',12);
ylabel('Overall SNR lost (linear)','FontSize',12);
grid minor;
% Modify colors
a = findobj(gca,'type','bar');
a(6).FaceColor = [0 104 255]./255;
a(1).FaceColor = [76 76 76]./255;
a(2).FaceColor = [127 127 127]./255;
a(3).FaceColor = [178 178 178]./255;
a(4).FaceColor = [204 204 204]./255;
a(5).FaceColor = [229 229 229]./255;
xlim([0.5 3.5]);
pos = get(gcf, 'Position');
set(gcf,'position',[pos(1),pos(2),677,235]);
grid minor;