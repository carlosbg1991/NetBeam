%% Multi-User Transmit Beamforming with USRP(R) Hardware
% This example uses beamforming techniques to send two different payloads
% to two receivers simultaneously in the same frequency band. Channel
% estimates from the receivers are used continuously to update the
% transmitted beams.
%
% The transmitter uses multiple USRP(R) radios as a composite radio with 4
% channels. Two X300/X310 radios or four N-series radios are required. All
% radios for the transmitter must be connected to the same PPS and 10 MHz
% clock generator via cables of equal lengths.
% 
% Besides the transmitter radios, the host computer must be connected to
% two more USRP(R) radios, one for each receiver. The receiver radios must
% also be connected to the same 10 MHz clock generator.
%
% To run this example, run MultiUserBeamformingExample, and then run
% helperMUBeamformRx1 and helperMUBeamformRx2 in two separate MATLAB
% sessions on the same host computer.
%
% Copyright 2016 The MathWorks, Inc.

%% Configure radios
clear all; close all; clc;
% % Find attached radios and allocate radios for transmitter and receivers.
% radioConfig = helperMUBeamformAllocateRadios;
% 
% % Save radio configuration to a file for helperMUBeamformRx1 and helperMUBeamformRx2
% save(fullfile(tempdir,'helperMUBeamformRadioConfig.mat'),'radioConfig');

% Configure experiment
numTxAntennas = 1;  % Select between 1, 2, 3 and 4
maxIter = 30000;

% Load Radio configurations
load(fullfile(tempdir,'helperMUBeamformRadioConfig.mat'),'radioConfig');

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
                                   'Gain',8,...
                                   'ClockSource','External',...
                                   'PPSSource','External');

%% Construct training signals and payloads
trainingSig = helperMUBeamformInitGoldSeq; % Based on Gold Sequences                    
trainingSig = trainingSig(:,1:numTxAntennas);  % Store only the ones being used

if ~exist(fullfile(tempdir,'BER_TxBits.mat'),'file')
    % Create sequence of bits for 64QAM
    M = 64;  k = log2(M);
    data{1} = randi([0 1],64*k,1);
    % Create sequence of bits for 32QAM
    M = 32;  k = log2(M);
    data{2} = randi([0 1],64*k,1);
    % Create sequence of bits for 16QAM
    M = 16;  k = log2(M);
    data{3} = randi([0 1],64*k,1);
    % Create sequence of bits for 8QAM
    M = 8;  k = log2(M);
    data{4} = randi([0 1],64*k,1);
    % Create sequence of bits for QPSK
    M = 4;  k = log2(M);
    data{5} = randi([0 1],64*k,1);
    % Create sequence of bits for BPSK
    M = 2;  k = log2(M);
    data{6} = randi([0 1],64*k,1);
    % Create sequence of bits for QPSK  - Extras
    M = 4;  k = log2(M);
    data_QPSK_extra = randi([0 1],64*k,2);
    % Save
    save(fullfile(tempdir,'BER_TxBits.mat'),'data','data_QPSK_extra');
else
    load(fullfile(tempdir,'BER_TxBits.mat'),'data','data_QPSK_extra');
end

% Create 64-QAM symbols
sym_64QAM = qammod(data{1},64,'InputType','bit','UnitAveragePower',true);
sym_64QAM = sym_64QAM.';
% Create 32-QAM symbols
sym_32QAM = qammod(data{2},32,'InputType','bit','UnitAveragePower',true);
sym_32QAM = sym_32QAM.';
% Create 16-QAM symbols
sym_16QAM = qammod(data{3},16,'InputType','bit','UnitAveragePower',true);
sym_16QAM = sym_16QAM.';
% Create 8-QAM symbols
sym_8QAM = qammod(data{4},8,'InputType','bit','UnitAveragePower',true);
sym_8QAM = sym_8QAM.';
% Create QPSK symbols
sym_QPSK = qammod(data{5},4,'InputType','bit','UnitAveragePower',true);
sym_QPSK = sym_QPSK.';
% Create QPSK symbols
sym_BPSK = qammod(data{6},2,'InputType','bit','UnitAveragePower',true);
sym_BPSK = sym_BPSK.';
% Create QPSK symbols
sym_QPSK_extra = qammod(data_QPSK_extra,4,'InputType','bit','UnitAveragePower',true);
sym_QPSK_extra = sym_QPSK_extra.';

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

%% Initialize variables and temporary files for channel feedback

% Have no knowledge of the channel yet - Generate channel estimate randomly
lastFeedback1 = rand(1,numTxAntennas) + 1j*rand(1,numTxAntennas);
fid1 = fopen(fullfile(tempdir,'helperMUBeamformfeedback1.bin'),'wb');
fwrite(fid1,[real(lastFeedback1) imag(lastFeedback1)],'double');
fclose(fid1);

% Open files for reading only
fid1 = fopen(fullfile(tempdir,'helperMUBeamformfeedback1.bin'),'rb');

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
        fprintf('h = %.7f + %.7fj\t\n',real(channelEst1(1)),imag(channelEst1(1)));
    end
    fprintf('\n');

    % Send signals to the radios
    txSig = [trainingSig; zeros(400,numTxAntennas); payload; zeros(100,numTxAntennas)] * 0.2;
    transmitter(txSig);
end

release(transmitter);

%% Appendix
% This example uses the following helper functions:
%
% * <matlab:edit('helperMUBeamformAllocateRadios.m') helperMUBeamformAllocateRadios.m>
% * <matlab:edit('helperMUBeamformInitGoldSeq.m') helperMUBeamformInitGoldSeq.m>
% * <matlab:edit('helperMUBeamformRx1.m') helperMUBeamformRx1.m>
% * <matlab:edit('helperMUBeamformRx2.m') helperMUBeamformRx2.m>
% * <matlab:edit('helperMUBeamformEstimateChannel.m') helperMUBeamformEstimateChannel.m>

displayEndOfDemoMessage(mfilename) 

%% Copyright Notice
% Universal Software Radio Peripheral(R) and USRP(R) are trademarks of
% National Instruments Corp.