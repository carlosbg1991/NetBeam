% % Configure radio
% radioConfig.Platform             = 'B210';
% radioConfig.SerialNum            = '30BC5C4';
% radioConfig.MasterClockRate      = 20e6;
% radioConfig.Decimationfactor     = 50;
% radioConfig.CenterFrequency      = 900e6;
% radioConfig.Gain                 = 8;
% radioConfig.SamplesPerFrame      = 200e3;
% radioConfig.EnableBurstMode      = true;
% radioConfig.NumFramesInBurst     = 1;
% radioConfig.OutputDataType       = 'double';
% radioConfig.ClockSource          = 'External'; % Synchronize transmitter and receiver in frequency
% % Connect to radio
% receiver = comm.SDRuReceiver('Platform',radioConfig.Platform,...
%                              'SerialNum',radioConfig.SerialNum,...
%                              'MasterClockRate',radioConfig.MasterClockRate,...
%                              'DecimationFactor',radioConfig.Decimationfactor,...
%                              'CenterFrequency',radioConfig.CenterFrequency,...
%                              'Gain',radioConfig.Gain,...
%                              'SamplesPerFrame',radioConfig.SamplesPerFrame,...
%                              'EnableBurstMode',radioConfig.EnableBurstMode,...
%                              'NumFramesInBurst',radioConfig.NumFramesInBurst,...
%                              'OutputDataType',radioConfig.OutputDataType,...
%                              'ClockSource',radioConfig.ClockSource);
% % Construct training signals - Golay Sequence 1 (do not modify)
% goldSeq = comm.GoldSequence;
% goldSeq.FirstPolynomial = [11 2 0];
% goldSeq.SecondPolynomial = [11 8 5 2 0];
% goldSeq.SamplesPerFrame = 2047;
% goldSeq.Index = 15;
% goldSeq.FirstInitialConditions = [0 1 0 1 1 0 1 0 0 1 1];
% goldSeq.SecondInitialConditions = [0 0 0 1 1 0 1 0 1 1 0];
% goldSequence1 = goldSeq();
% txfilter = comm.RaisedCosineTransmitFilter;
% goldSeqRef = txfilter([goldSequence1 ; zeros(10,1)]); % Pad extra zeros at the end
% % Configuration
% L = size(goldSeqRef,1);
% % Main loop
% maxIter = 200;  % Maximum number of iterations
% % SNR     = 15;  % Signal-To-Noise Ratio
% chRealTot = zeros(maxIter,1);
% chEstiTot = zeros(maxIter,1);
for iter = 1:maxIter
    [rxSig, len] = receiver();
%     chRealTot(iter)  = rand(1,1) + 1i.*rand(1,1);
%     rxSig = [zeros(5000,1);chRealTot(iter).*goldSeqRef;zeros(5000,1)];
%     rxPower = sum(var(chRealTot(iter).*goldSeqRef));
%     SNRlin = db2pow(SNR);
%     noise = rand(length(rxSig),1) + 1i.*rand(length(rxSig),1);
%     noisePower = sum(var(noise));
%     noise = sqrt(rxPower/(SNRlin*noisePower)).*noise;
%     rxSig = rxSig + noise;  % Add noise to signal
%     len = length(rxSig);
    if len > 0
        % Estimate the channel
        rxSigLen = length(rxSig);
        crossCorr = xcorr(rxSig,goldSeqRef(:,1));
        crossCorr = crossCorr((rxSigLen+1):end);
        peakIntervals = find(abs(crossCorr)>(0.8*max(abs(crossCorr)))); % At least 80% of global maximum
        % We know that the training signal has more than 1000 samples
        peakCandidates = find((peakIntervals(2:end)-peakIntervals(1:end-1))>1000);
%         peakCandidates = (1:1:length(peakIntervals));
        % Don't use the last peak found
        if ~isempty(peakCandidates);  peakCandidates(end) = [];  end
        channelEstimate = zeros(length(peakCandidates),1);
        for idxPeak = 1:length(peakCandidates)
            startIndex = peakIntervals(peakCandidates(idxPeak)+1);
            % Focus on a small window that contains a peak. +/- 100 samples
            if startIndex > 100;  windowRange = startIndex+(-100:100);
            else;                 windowRange = 1:(startIndex+100);
            end
            dataWindow = abs(crossCorr(windowRange));
            % Find exact locations of the peaks (one in each of 4 cross correlations)
            [~,maxLoc] = max(dataWindow); % maxLoc is a 1x4 vector
            finalPeakLoc = startIndex-101;
            indexAll = finalPeakLoc + maxLoc + rxSigLen;
            indexMin = min(indexAll);
            indexMax = max(indexAll);
            % Retrieve received training signal
            receivedTrainingSig = rxSig((indexMin-rxSigLen+1):(indexMax-rxSigLen+L));
            % Use least squares fit to estimate channel response
            channelEstimate(idxPeak) = goldSeqRef\receivedTrainingSig;
        end
    end
    chEstiTot(iter) = mean(channelEstimate.');
%     fprintf('iter %d\n',iter);
%     fprintf('\tReal channel: %.3f + %.3fj\n',real(chRealTot(iter)),imag(chRealTot(iter)));
	fprintf('\titer %d - h_est: %.7f + %.7fj\n',iter,real(chEstiTot(iter)),imag(chEstiTot(iter)));
end

figure;
subplot(311); hold on; grid minor;
% plot((1:maxIter),real(chRealTot),'Color','b');
plot((1:maxIter),real(chEstiTot),'Color','r','LineWidth',2);
title('Real','FontSize',12);
xlabel('Iteration','FontSize',12);
ylabel('Gain','FontSize',12);
subplot(312); hold on; grid minor;
% plot((1:maxIter),imag(chRealTot),'Color','b');
plot((1:maxIter),imag(chEstiTot),'Color','r','LineWidth',2);
title('Imaginary','FontSize',12);
xlabel('Iteration','FontSize',12);
ylabel('Gain','FontSize',12);
subplot(313); hold on; grid minor;
plot((1:maxIter),abs(chEstiTot),'Color','r','LineWidth',2);
title('Absolute Gain','FontSize',12);
xlabel('Iteration','FontSize',12);
ylabel('Gain','FontSize',12);

figure; hold on;
covTot_imag = zeros(200,1);
covTot_real = zeros(200,1);
re = real(chEstiTot);
im = imag(chEstiTot);
re(isnan(re))=0;
window_size = 1;
for k = 1:200-window_size
    covTot_imag(k) = cov(im(k:k+window_size));
    covTot_real(k) = cov(re(k:k+window_size));
end
plot(abs(1-covTot_imag),'LineWidth',2);
plot(abs(1-covTot_real),'LineWidth',2);
legend('Imaginary','Real');

release(receiver);

save('sim_TxStat_RxMob.mat');

%% Appendix
% This example uses the following helper functions:
%
% * <matlab:edit('helperMUBeamformInitGoldSeq.m') helperMUBeamformInitGoldSeq.m>
% * <matlab:edit('helperMUBeamformEstimateChannel.m') helperMUBeamformEstimateChannel.m>

%% Copyright Notice
% Universal Software Radio Peripheral(R) and USRP(R) are trademarks of
% National Instruments Corp.