%% Multi-User Transmit Beamforming with USRP(R) Hardware
% Companion script for MultiUserBeamformingExample. Run this script after
% running MultiUserBeamformingExample in a separate MATLAB session.
% See details in MultiUserBeamformingExample.m
%
% Copyright 2016 The MathWorks, Inc.

% % Load radio configuration from a file
% if exist(fullfile(tempdir,'helperMUBeamformRadioConfig.mat'),'file') == 0
%     error('MultiUserBeamformingExample must be running in a separate MATLAB session on this computer first.');
% end
% load(fullfile(tempdir,'helperMUBeamformRadioConfig.mat'));
% 
% % Connect to radio
% receiver = comm.SDRuReceiver('Platform',radioConfig.rxPlatform1,...
%                              radioConfig.rxIDProp1,radioConfig.rxID1);
% receiver.MasterClockRate = radioConfig.rxMasterClockRate1;
% receiver.DecimationFactor = radioConfig.rxDecimationfactor1;
% receiver.ClockSource = 'External'; % Synchronize transmitter and receiver in frequency
% receiver.Gain = 8;
% receiver.CenterFrequency = 900e6;
% receiver.SamplesPerFrame = 200e3;
% receiver.EnableBurstMode = true;
% receiver.NumFramesInBurst = 1;
% receiver.OutputDataType = 'double';
% 
% % Radio settings
% receiver
% 
% % Open temporary file for writing channel estimate
% fid = fopen(fullfile(tempdir,'helperMUBeamformfeedback1.bin'),'wb');
% 
% % Get Gold sequences for estimating channel response
% goldSeqRef = helperMUBeamformInitGoldSeq;
% 
% % Variable to store the BER
% BER = zeros(5000,6);
% 
% % Load Transmitted bits for the modulations used
% load(fullfile(tempdir,'BER_TxBits.mat'),'data');
% 
% tic;
% maxIter = 5000;
% modList = [64 32 16 8 4 2];
% chTot = zeros(4,maxIter);
% % Main loop
for i = 1:maxIter
    elapsedTime = toc;
    [rxSig, len] = receiver();
    if len > 0
        [channelEstimate, payload] = ...
            helperMUBeamformEstimateChannel(rxSig, goldSeqRef);
        fftOut = fft(reshape(payload, 256, 64));
        
        for modIdx = 1:length(modList)
            index = 4 + modIdx;  % First 4 subcarriers contain 0's
            y = fftOut(index,:).';  % Extract Subcarrier
            y = y/sqrt(mean(y'*y));
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
        
        chTot(:,i) = channelEstimate;
    end
    elapsedOld = elapsedTime;
    elapsedTime = toc;
%     fprintf('Total Elapsed:  %.3f\n',elapsedTime);
%     fprintf('Iteration time: %.3f\n',elapsedTime - elapsedOld);
    fprintf('\titer %d - h_est(1): %.7f + %.7fj\n',i,real(chTot(1,i)),imag(chTot(1,i)));
end

save('sim_BER-MU_TxAll_RxAll.mat');

lastIter = i;
for idx = 1:4
    figure(15);
    subplot(3,4,idx); hold on; grid minor;
    % plot((1:maxIter),real(chRealTot),'Color','b');
    plot((1:lastIter),real(chTot(idx,1:lastIter)),'Color','r','LineWidth',2);
    title('Real','FontSize',12);
    xlabel('Iteration','FontSize',12);
    ylabel('Gain','FontSize',12);
    subplot(3,4,idx+4); hold on; grid minor;
    % plot((1:maxIter),imag(chRealTot),'Color','b');
    plot((1:lastIter),imag(chTot(idx,1:lastIter)),'Color','r','LineWidth',2);
    title('Imaginary','FontSize',12);
    xlabel('Iteration','FontSize',12);
    ylabel('Gain','FontSize',12);
    subplot(3,4,idx+8); hold on; grid minor;
    plot((1:lastIter),abs(chTot(idx,1:lastIter)),'Color','r','LineWidth',2);
    title('Absolute Gain','FontSize',12);
    xlabel('Iteration','FontSize',12);
    ylabel('Gain','FontSize',12);

    figure(16);
    subplot(4,1,idx); hold on;
    covTot_imag = zeros(lastIter,1);
    covTot_real = zeros(lastIter,1);
    re = real(chTot(idx,:));
    im = imag(chTot(idx,:));
    re(isnan(re))=0;
    window_size = 1;
    for k = 1:1:lastIter-window_size
        covTot_imag(k) = cov(im(k:k+window_size));
        covTot_real(k) = cov(re(k:k+window_size));
    end
    plot(abs(1-covTot_imag),'LineWidth',2);
    plot(abs(1-covTot_real),'LineWidth',2);
    legend('Imaginary','Real');
    ylim([0.999 1]);
end

figure(17);
subplot(611); plot(abs(BER(1:i,1))); legend('64-QAM'); ylabel('BER'); ylim([0 1]);
subplot(612); plot(abs(BER(1:i,2))); legend('32-QAM'); ylabel('BER'); ylim([0 1]);
subplot(613); plot(abs(BER(1:i,3))); legend('16-QAM'); ylabel('BER'); ylim([0 1]);
subplot(614); plot(abs(BER(1:i,4))); legend('8-QAM'); ylabel('BER'); ylim([0 1]);
subplot(615); plot(abs(BER(1:i,5))); legend('QPSK'); ylabel('BER'); ylim([0 1]);
subplot(616); plot(abs(BER(1:i,6))); legend('BPSK'); ylabel('BER'); ylim([0 1]);

release(receiver);

%% Appendix
% This example uses the following helper functions:
%
% * <matlab:edit('helperMUBeamformInitGoldSeq.m') helperMUBeamformInitGoldSeq.m>
% * <matlab:edit('helperMUBeamformEstimateChannel.m') helperMUBeamformEstimateChannel.m>

%% Copyright Notice
% Universal Software Radio Peripheral(R) and USRP(R) are trademarks of
% National Instruments Corp.