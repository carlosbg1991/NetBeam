%% Multi-User Transmit Beamforming with USRP(R) Hardware

% Configuration
maxIter = 110;
modList = [64 32 16 8 4 2];
nTxAntennas = 1;

% % Check for required files
% if exist(fullfile(tempdir,'helperMUBeamformRadioConfig.mat'),'file') == 0
%     error('MultiUserBeamformingExample must be running in a separate MATLAB session on this computer first.');
% elseif exist(fullfile(tempdir,'helperMUBeamformRadioConfig.mat'),'file') == 0
%     error('Need to run BER_MultiUserBeamformingExample first to generate tx symbols.');
% end
% % Load radio configuration from a file
% load(fullfile(tempdir,'helperMUBeamformRadioConfig.mat'));
% % Load Transmitted bits for the modulations used
% load(fullfile(tempdir,'BER_TxBits.mat'),'data');
% 
% % Connect to radio
% receiver = comm.SDRuReceiver('Platform',radioConfig.rxPlatform1,...
%                              radioConfig.rxIDProp1,radioConfig.rxID1,...
%                              'MasterClockRate',radioConfig.rxMasterClockRate1,...
%                              'DecimationFactor',radioConfig.rxDecimationfactor1,...
%                              'ClockSource','External',...
%                              'Gain',8,...
%                              'CenterFrequency',900e6,...
%                              'EnableBurstMode',true,...
%                              'SamplesPerFrame',200e3,...
%                              'NumFramesInBurst',1,...
%                              'OutputDataType','double');
%                          
% % Open temporary file for writing channel estimate
% fid = fopen(fullfile(tempdir,'helperMUBeamformfeedback1.bin'),'wb');
% 
% % Get Gold sequences for estimating channel response
% goldSeqRef = helperMUBeamformInitGoldSeq;
% 
% % Variable to store the BER
% BER = zeros(5000,6);

tic;
chTot = zeros(nTxAntennas,maxIter);
% Main loop
for i = 1:maxIter
    elapsedTime = toc;
    [rxSig, len] = receiver();
    if len > 0
        [channelEstimate, payload] = ...
            BER_helperMUBeamformEstimateChannel(rxSig, goldSeqRef, nTxAntennas);
        fftOut = fft(reshape(payload, 256, 64));
        
        for modIdx = 1:length(modList)
            index = 4 + modIdx;  % First 4 subcarriers contain 0's
            y = fftOut(index,:).';  % Extract Subcarrier
            y = y/sqrt(mean(y'*y));
            y = 1/sqrt(sum(var(y))).*y;  % Normalize symbols
            % Compensate for impairments
%             [y,~] = freqComp(y);  % Compensate for Frequency offset
%             y = iqImbComp(y);  % Compensate for IQ Imbalance
            % Plot constellation
            figure(modIdx); clf('reset');
            figure(modIdx); hold on;
            y_tx = qammod(data{modIdx},modList(modIdx),'InputType','bit','UnitAveragePower',true);
            plot(real(y_tx),imag(y_tx),'LineStyle','None','Marker','.','Color','r');
            plot(real(y),imag(y),'LineStyle','None','Marker','.','Color','b');
            xlim([-2 2]);  ylim([-2 2]);  % Normalized
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
    fprintf('Iter %d:\t',i);
    for id = 1:nTxAntennas
        fprintf('h = %.7f + %.7fj\t',real(chTot(1,i)),imag(chTot(1,i)));
    end
    fprintf('\n');
end

save('sim_BER-exp2ant.mat');

lastIter = i;
for idx = 1:nTxAntennas
    figure(15);
    subplot(nTxAntennas,3,3*(idx-1) + 1); hold on; grid minor;
    % plot((1:maxIter),real(chRealTot),'Color','b');
    plot((1:lastIter),real(chTot(idx,1:lastIter)),'Color','r','LineWidth',2);
    title('Real','FontSize',12);
    xlabel('Iteration','FontSize',12);
    ylabel('Gain','FontSize',12);
    subplot(nTxAntennas,3,3*(idx-1) + 2); hold on; grid minor;
    % plot((1:maxIter),imag(chRealTot),'Color','b');
    plot((1:lastIter),imag(chTot(idx,1:lastIter)),'Color','r','LineWidth',2);
    title('Imaginary','FontSize',12);
    xlabel('Iteration','FontSize',12);
    ylabel('Gain','FontSize',12);
    subplot(nTxAntennas,3,3*(idx-1) + 3); hold on; grid minor;
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
totGain = 0;
for txID = 1:nTxAntennas
    t = abs(chTot(txID,1:lastIter));
    t(isnan(t))=0;  % replace nan values
    totGain = totGain + t;
end
 plot((1:lastIter),totGain,'Color','r','LineWidth',2);
title('Total channel Gain','FontSize',12);
xlabel('Iteration','FontSize',12);
ylabel('Gain (linear)','FontSize',12);

figure(18);
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