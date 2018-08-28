%% Clear environment
clear all; close all; clc;

%% Configure experiment
fileName = 'data/channelEstimation.bin';  % File location for channel estimation
numTxAntennas = 4;  % Select between 1, 2 and 4
maxIter = 30000;  % Maximum transmissions over the air
gain = -10; % in dB
NFFT = 256;  % Point for the FFT (OFDM modulator)
Ndc = 6;  % Zero padding for central frequencies
Nedge = 10;  % Zero padding for the edges (multipath effect)

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
    bits{1} = randi([0 1],64*k,1);
    M = 32;  k = log2(M);  % 32QAM
    bits{2} = randi([0 1],64*k,1);
    M = 16;  k = log2(M);  % 16QAM
    bits{3} = randi([0 1],64*k,1);
    M = 8;  k = log2(M);  % 8QAM
    bits{4} = randi([0 1],64*k,1);
    M = 4;  k = log2(M);  % QPSK
    bits{5} = randi([0 1],64*k,1);
    M = 2;  k = log2(M);  % BPSK
    bits{6} = randi([0 1],64*k,1);
    
    symbols(1,:) = qammod(bits{1},64,'InputType','bit','UnitAveragePower',true).';  % 64QAM
    symbols(2,:) = qammod(bits{2},32,'InputType','bit','UnitAveragePower',true).';  % 32QAM
    symbols(3,:) = qammod(bits{3},16,'InputType','bit','UnitAveragePower',true).';  % 16QAM
    symbols(4,:)= qammod(bits{4},8,'InputType','bit','UnitAveragePower',true).';  % 8QAM
    symbols(5,:)= qammod(bits{5},4,'InputType','bit','UnitAveragePower',true).';  % QPSK
    symbols(6,:)= qammod(bits{6},2,'InputType','bit','UnitAveragePower',true).';  % BPSK
    
    %% Construct payload 1:
    % 64 symbols with IFFT length of 256
    % Each symbol uses 8 subcarriers
    % Subcarrier 5 uses 64-QAM. Each point of 64-QAM is used once.
    % Other subcarriers use QPSK
    BlockSize = (NFFT/2) - (Ndc/2) - Nedge;
    rep = floor(BlockSize/6);
    modOut1 = zeros(NFFT,64);
    idxFFT(:,1) = [(Ndc/2)+0*rep+1:(Ndc/2)+1*rep (NFFT/2)+0*rep+1:(NFFT/2)+1*rep];  % 64QAM
    idxFFT(:,2) = [(Ndc/2)+1*rep+1:(Ndc/2)+2*rep (NFFT/2)+1*rep+1:(NFFT/2)+2*rep];  % 32QAM
    idxFFT(:,3) = [(Ndc/2)+2*rep+1:(Ndc/2)+3*rep (NFFT/2)+2*rep+1:(NFFT/2)+3*rep];  % 16QAM
    idxFFT(:,4) = [(Ndc/2)+3*rep+1:(Ndc/2)+4*rep (NFFT/2)+3*rep+1:(NFFT/2)+4*rep];  % 8QAM
    idxFFT(:,5) = [(Ndc/2)+4*rep+1:(Ndc/2)+5*rep (NFFT/2)+4*rep+1:(NFFT/2)+5*rep];  % QPSK
    idxFFT(:,6) = [(Ndc/2)+5*rep+1:(Ndc/2)+6*rep (NFFT/2)+5*rep+1:(NFFT/2)+6*rep];  % BPSK
    modOut1(idxFFT(:,1),:) = repmat(symbols(1,:),2*rep,1);
    modOut1(idxFFT(:,2),:) = repmat(symbols(2,:),2*rep,1);
    modOut1(idxFFT(:,3),:) = repmat(symbols(3,:),2*rep,1);
    modOut1(idxFFT(:,4),:)  = repmat(symbols(4,:),2*rep,1);
    modOut1(idxFFT(:,5),:)  = repmat(symbols(5,:),2*rep,1);
    modOut1(idxFFT(:,6),:)  = repmat(symbols(6,:),2*rep,1);
    % OFDM Modulation and append OFDM blocks to create time frame
    payload = reshape(ifft(modOut1),[],1);
    % Scale time-domain signal appropriately
    payload = payload/max(real(payload))*0.5;

    % Save
    bits{1} = repmat(bits{1},2*rep,1);
    bits{2} = repmat(bits{1},2*rep,1);
    bits{3} = repmat(bits{1},2*rep,1);
    bits{4} = repmat(bits{1},2*rep,1);
    bits{5} = repmat(bits{1},2*rep,1);
    bits{6} = repmat(bits{1},2*rep,1);
    save(fullfile('data','information.mat'),'bits','symbols','idxFFT');
else
    load(fullfile('data','information.mat'),'bits','symbols','idxFFT');
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
    payload = payload*beamWeight1;
    
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