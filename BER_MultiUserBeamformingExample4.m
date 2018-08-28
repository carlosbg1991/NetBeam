%% Clear environment
clear all; close all; clc;

%% Configure experiment
fileName = 'data/channelEstimation.bin';  % File location for channel estimation
numTxAntennas = 4;  % Select between 1, 2 and 4
maxIter = 30000;  % Maximum transmissions over the air
gain = -10; % in dB

%% Configure radios

% Load Radio configurations
load(fullfile('data','radioConfig.mat'));

% Configure transmitters
txIPAvail = split(radioConfig.txIPAddrs,',');
txIPAvail = {txIPAvail{2},txIPAvail{1}}.';
switch numTxAntennas
    case 1
        radioConfig.ChannelMapping = 1;  % Use only 1 antenna
        radioConfig.txIPAddrs = txIPAvail{1};  % only 1 radio
    case 2
        radioConfig.ChannelMapping = [1 2];  % Use 2 antennas from the same USRP
        radioConfig.txIPAddrs = txIPAvail{1};  % only 1 radio
    case 3
        radioConfig.ChannelMapping = [1 2 3];  % Use 3 antennas
	case 4
        radioConfig.ChannelMapping = [1 2 3 4];  % Use 4 antennas
    otherwise
        error('Wrong number of tx antennas for numTxAntennas');
end

% Create Transmitter SO
transmitter = comm.SDRuTransmitter('Platform',radioConfig.txPlatform,...
                                   'IPAddress',radioConfig.txIPAddrs,...
                                   'MasterClockRate',radioConfig.txMasterClockRate,...
                                   'InterpolationFactor',radioConfig.txInterpolationfactor,...
                                   'ChannelMapping',radioConfig.ChannelMapping,...
                                   'CenterFrequency',900e6,...
                                   'Gain',gain,...
                                   'ClockSource','External',...
                                   'PPSSource','External');

%% Construct training signals and payloads
if ~exist(fullfile('data','trainingSig.mat'),'file')
    % Generate Gold sequences as training signals
    trainingSig = helperMUBeamformInitGoldSeq;
    % Save pre-defined Gold sequences
    save(fullfile('data','trainingSig.mat'),'trainingSig');
else
    % Load pre-defined Gold sequences
    load(fullfile('data','trainingSig.mat'),'trainingSig');
end
trainingSig = trainingSig(:,1:numTxAntennas);  % Store only the ones being used

%% Create data
if ~exist(fullfile('data','information.mat'),'file')
    M = 64;  k = log2(M);  % 64QAM
    data{1} = randi([0 1],64*k,1);
    M = 32;  k = log2(M);  % 32QAM
    data{2} = randi([0 1],64*k,1);
    M = 16;  k = log2(M);  % 16QAM
    data{3} = randi([0 1],64*k,1);
    M = 8;  k = log2(M);  % 8QAM
    data{4} = randi([0 1],64*k,1);
    M = 4;  k = log2(M);  % QPSK
    data{5} = randi([0 1],64*k,1);
    M = 2;  k = log2(M);  % BPSK
    data{6} = randi([0 1],64*k,1);
    M = 4;  k = log2(M);  % BPSK - extra
    data_QPSK_extra = randi([0 1],64*k,2);
    
    % Create 64-QAM symbols
    sym_64QAM = qammod(data{1},64,'InputType','bit','UnitAveragePower',true).';
    % Create 32-QAM symbols
    sym_32QAM = qammod(data{2},32,'InputType','bit','UnitAveragePower',true).';
    % Create 16-QAM symbols
    sym_16QAM = qammod(data{3},16,'InputType','bit','UnitAveragePower',true).';
    % Create 8-QAM symbols
    sym_8QAM = qammod(data{4},8,'InputType','bit','UnitAveragePower',true).';
    % Create QPSK symbols
    sym_QPSK = qammod(data{5},4,'InputType','bit','UnitAveragePower',true).';
    % Create QPSK symbols
    sym_BPSK = qammod(data{6},2,'InputType','bit','UnitAveragePower',true).';
    % Create QPSK symbols
    sym_QPSK_extra = qammod(data_QPSK_extra,4,'InputType','bit','UnitAveragePower',true).';

    % Construct payload 1:
    % 64 symbols with IFFT length of 256
    % Each symbol uses 8 subcarriers
    % Subcarrier 5 uses 64-QAM. Each point of 64-QAM is used once.
    % Other subcarriers use QPSK
    modOut1 = [zeros(4,64);       % Zero padding
               sym_64QAM;         % 64-QAM - 1 subcarrier
               sym_32QAM;         % 32-QAM - 1 subcarrier
               sym_16QAM;         % 16-QAM - 1 subcarrier
               sym_8QAM;          % 8-QAM - 1 subcarrier
               sym_QPSK;          % QPSK - 1 subcarrier
               sym_BPSK;          % BPSK - 1 subcarrier
               sym_QPSK_extra;    % QPSK - 1 subcarrier - Extra
               zeros(256-4-8,64)]; % Zero padding
    payload1 = reshape(ifft(modOut1),[],1);
    % Scale time-domain signal appropriately
    payload1 = payload1/max(real(payload1))*0.5;
    % Save
    save(fullfile('data','information.mat'),'data','data_QPSK_extra');
else
    load(fullfile('data','information.mat'),'data','data_QPSK_extra');
end

%% Initialize variables and temporary files for channel feedback

% Initialize channel estimate randomly - create file
fid = fopen(fileName,'wb');
chEst_prelim = rand(1,numTxAntennas) + 1j*rand(1,numTxAntennas);
fwrite(fid,[real(chEst_prelim) imag(chEst_prelim)],'double');
fclose(fid);

% Open file to share channel estimation(read only)
fid1 = fopen(fileName,'rb');

%% Main loop

disp('Sending ONE payload to ONE receivers simultaneously in the same frequency band ...')
disp('Channel estimates from the receiver are used continuously to update the transmitted beams.')
disp('Please run helperMUBeamformRx1_3 in A SEPARATE MATLAB sessions on this computer.')

for i = 1:maxIter
    fseek(fid1,0,'bof'); % Read from the beginning of the file
    % The channel between 4 TX antennas and 1 RX antenna is modeled
    % by 4 complex gains. This approximation works because the
    % signal has very narrow bandwidth (400k samples per second).
    channelEst1 = fread(fid1,numTxAntennas*2,'double');
    if length(channelEst1) == numTxAntennas*2
        % File content has expected length
        channelEst1 = (   channelEst1(1:numTxAntennas) + ...
                       1j*channelEst1(numTxAntennas+1:numTxAntennas*2)).';
        lastFeedback1 = channelEst1;
    else
        % Use last feedback
        channelEst1 = lastFeedback1;
    end

    % Rotate payload 1 so that receiver 1 does not need to
    % correct the phase of payload 1
    beamWeight1 = ones(1,numTxAntennas);  % Does not matter
    phaseCorrection1 = beamWeight1 * channelEst1.';  % Beamforming info
    phaseCorrection1 = phaseCorrection1/abs(phaseCorrection1);
    beamWeight1 = beamWeight1 ./ phaseCorrection1;

    % Beamforming for payload
    payload = payload1*beamWeight1;
    
    fprintf('Iter %d - Applying channel:\n',i);
    for id = 1:numTxAntennas
        fprintf('h = %.7f + %.7fj\t\n',real(channelEst1(id)),imag(channelEst1(id)));
    end
    fprintf('\n');

    % Send signals to the radios
    txSig = [trainingSig; zeros(400,numTxAntennas); payload; zeros(100,numTxAntennas)] * 0.2;
    transmitter(txSig);
end

release(transmitter);  % Release the System Object for future use