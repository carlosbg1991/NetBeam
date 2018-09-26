function BER_MultiUserBeamformingExample4(varargin)
% Multi-User Beamforming transmitter with the following features:
% - Local file to store the channel estimation.
% - Only 6 subcarriers are used for data transmision.
% - The number of antennas and gain are variable.

%% Configure experiment
if (nargin==4)
    numTxAntennas = varargin{1};
    maxIter       = varargin{2};
    gain          = varargin{3};
    fileName      = varargin{4};
elseif (nargin==0)
    numTxAntennas = 4;  % Select between 1, 2 and 4
    maxIter       = 30000;  % Maximum transmissions over the air
    gain          = 0; % in dB
    fileName      = 'mierda.bin';  % File location for channel estimation (decentralized)
%     fileName      = 'helperMUBeamformfeedback1.bin';  % File location for channel estimation (centralized)
else
    error('ERROR: The number of input arguments mismatch the expected.\n');
end

%% Configure radios

% Load Radio configurations
load(fullfile('data','radioConfig.mat'),'radioConfig');

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
if ~exist(fullfile('data','information4.mat'),'file')
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
    
    % Create 64-QAM symbols
    symbols(1,:) = qammod(bits{1},64,'InputType','bit','UnitAveragePower',true).';  % 64QAM
    symbols(2,:) = qammod(bits{2},32,'InputType','bit','UnitAveragePower',true).';  % 32QAM
    symbols(3,:) = qammod(bits{3},16,'InputType','bit','UnitAveragePower',true).';  % 16QAM
    symbols(4,:) = qammod(bits{4},8,'InputType','bit','UnitAveragePower',true).';  % 8QAM
    symbols(5,:) = qammod(bits{5},4,'InputType','bit','UnitAveragePower',true).';  % QPSK
    symbols(6,:) = qammod(bits{6},2,'InputType','bit','UnitAveragePower',true).';  % BPSK

    % Construct payload 1:
    % 64 symbols with IFFT length of 256
    % Each symbol uses 8 subcarriers
    % Subcarrier 5 uses 64-QAM. Each point of 64-QAM is used once.
    % Other subcarriers use QPSK
    modOut1 = [zeros(4,64);       % Zero padding
               symbols(1,:);         % 64-QAM - 1 subcarrier
               symbols(2,:);         % 32-QAM - 1 subcarrier
               symbols(3,:);         % 16-QAM - 1 subcarrier
               symbols(4,:);          % 8-QAM - 1 subcarrier
               symbols(5,:);          % QPSK - 1 subcarrier
               symbols(6,:);          % BPSK - 1 subcarrier
               zeros(256-4-6,64)]; % Zero padding
    payload1 = reshape(ifft(modOut1),[],1);
    % Scale time-domain signal appropriately
    payload1 = payload1/max(real(payload1))*0.5;
    % Save
    save(fullfile('data','information4.mat'),'bits','symbols','payload1');
else
    load(fullfile('data','information4.mat'),'bits','symbols','payload1');
end

%% Initialize variables and temporary files for channel feedback

% Open file to share channel estimation(read only)
fid1 = fopen(fileName,'rb');
lastFeedback1 = randn(1,numTxAntennas) + 1i.*randn(1,numTxAntennas);

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


% EOF