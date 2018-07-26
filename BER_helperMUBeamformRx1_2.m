%% Multi-User Transmit Beamforming with USRP(R) Hardware
% Companion script for MultiUserBeamformingExample. Run this script after
% running MultiUserBeamformingExample in a separate MATLAB session.
% See details in MultiUserBeamformingExample.m
%
% Copyright 2016 The MathWorks, Inc.

% Load radio configuration from a file
if exist(fullfile(tempdir,'helperMUBeamformRadioConfig.mat'),'file') == 0
    error('MultiUserBeamformingExample must be running in a separate MATLAB session on this computer first.');
end
load(fullfile(tempdir,'helperMUBeamformRadioConfig.mat'));

% Connect to radio
receiver = comm.SDRuReceiver('Platform',radioConfig.rxPlatform1,...
                             radioConfig.rxIDProp1,radioConfig.rxID1);
receiver.MasterClockRate = radioConfig.rxMasterClockRate1;
receiver.DecimationFactor = radioConfig.rxDecimationfactor1;
receiver.ClockSource = 'External'; % Synchronize transmitter and receiver in frequency
receiver.Gain = 8;
receiver.CenterFrequency = 900e6;
receiver.SamplesPerFrame = 200e3;
receiver.EnableBurstMode = true;
receiver.NumFramesInBurst = 1;
receiver.OutputDataType = 'double';

% Radio settings
receiver

% Set up spectrum analyzer and constellation diagram
spAnalyzer = dsp.SpectrumAnalyzer;
spAnalyzer.SampleRate = 400e3;
spAnalyzer.SpectralAverages = 64;
spAnalyzer.Title = 'Receiver 1';
spAnalyzer.YLimits = [-100 -20];

% Open temporary file for writing channel estimate
fid = fopen(fullfile(tempdir,'helperMUBeamformfeedback1.bin'),'wb');

% Get Gold sequences for estimating channel response
goldSeqRef = helperMUBeamformInitGoldSeq;

% Variable to store the BER
BER = zeros(5000,6);

% Load Transmitted bits for the modulations used
load(fullfile(tempdir,'BER_TxBits.mat'),'data');

tic;
modList = [64 32 16 8 4 2];
% Main loop
for i = 1:5000
    elapsedTime = toc;
    [rxSig, len] = receiver();
    if len > 0
        [channelEstimate, payload] = ...
            helperMUBeamformEstimateChannel(rxSig, goldSeqRef);
        spAnalyzer(payload);    % Plot power spectrum
        fftOut = fft(reshape(payload, 256, 64));
        
        for modIdx = 1:length(modList)
            index = 4 + modIdx;  % First 4 subcarriers contain 0's
            y = fftOut(index,:).';  % Extract Subcarrier
            y = y/sqrt(mean(y'*y));
%             constDiagram(y);        % Plot constellation
            y = 1/sqrt(sum(var(y))).*y;  % Normalize symbols
            % Plot constellation
            figure(modIdx); clf('reset');
            figure(modIdx); hold on;
            y_tx = qammod(data{modIdx},modList(modIdx),'InputType','bit','UnitAveragePower',true);
            plot(real(y_tx),imag(y_tx),'LineStyle','None','Marker','.','Color','r');
            plot(real(y),imag(y),'LineStyle','None','Marker','.','Color','b');
            tit = strcat('Receiver with k =',{' '},num2str(modList(modIdx)));
            title(tit{1},'FontSize',12);
            % Compute Bit Error Rate for the 64-QAM modulation
            if ~isempty(y) && ~any(isnan(y))
                % Demodulator expecting normalized symbols
                M = modList(modIdx);
                data_rx = qamdemod(y,M,'OutputType','bit','UnitAveragePower',true);
                BER(i,modIdx) = sum(abs(data{modIdx} - data_rx))/length(data_rx);
                fprintf('Iter %d - BER: %.3f\n',i,BER(i,modIdx));
            else
                BER(i,modIdx) = BER(i-1,modIdx);  % First element is the BER
                fprintf('Iter %d - BER: %.3f (Hardcoded)\n',i,BER(i,modIdx));
            end
        end
        
        % Write channel estimate to a file
        fseek(fid,0,'bof');
        if ~any(isnan(channelEstimate))
            fwrite(fid,[real(channelEstimate) imag(channelEstimate)],'double');
        end
    end
    elapsedOld = elapsedTime;
    elapsedTime = toc;
    fprintf('Total Elapsed:  %.3f\n',elapsedTime);
    fprintf('Iteration time: %.3f\n',elapsedTime - elapsedOld);
end

release(receiver);

%% Appendix
% This example uses the following helper functions:
%
% * <matlab:edit('helperMUBeamformInitGoldSeq.m') helperMUBeamformInitGoldSeq.m>
% * <matlab:edit('helperMUBeamformEstimateChannel.m') helperMUBeamformEstimateChannel.m>

%% Copyright Notice
% Universal Software Radio Peripheral(R) and USRP(R) are trademarks of
% National Instruments Corp.