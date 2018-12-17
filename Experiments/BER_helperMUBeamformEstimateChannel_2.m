function [channelEstimate, payload, finalCorrectionTime] = BER_helperMUBeamformEstimateChannel_2(rxSig, goldSeqRef, varargin)

if nargin == 2
    % Call regular channel estimator
    [channelEstimate, payload] = helperMUBeamformEstimateChannel(rxSig, goldSeqRef);
    return;
elseif nargin~=4
    error('Need to state the number of antennas used in Tx1 and Tx2');
end

nTxAntennas1 = varargin{1};
nTxAntennas2 = varargin{2};
nTxAntennas = nTxAntennas1 + nTxAntennas2;
rxSigLen = length(rxSig);

figure(1); cla reset; hold on;

% Detect training signal in TX1 by cross correlation
crossCorr = cell(nTxAntennas,1);
for txID = 1:nTxAntennas1
    crossCorr{txID} = xcorr(rxSig,goldSeqRef(:,txID));
    crossCorr{txID} = crossCorr{txID}((rxSigLen+1):end);
    plot(abs(crossCorr{txID}));
end

% Detect training signal in TX2 by cross correlation
for txID = 1:nTxAntennas2
    crossCorr{nTxAntennas1+txID} = xcorr(rxSig,goldSeqRef(:,2+txID));
    crossCorr{nTxAntennas1+txID} = crossCorr{nTxAntennas1+txID}((rxSigLen+1):end);
    plot(abs(crossCorr{nTxAntennas1+txID}));
end

L = size(goldSeqRef,1);
payloadLen = 64*256; % 64 symbols with IFFT length of 256
payload = zeros(payloadLen,1);
numPayload = 0;
finalCorrectionTime = [0 0];  % Initialize the correction time
myStartIdx1 = [];

%% FOR TX1
if nTxAntennas1 > 0
    % Find intervals containing peaks of cross correlation FOR CHANNEL 1 - TX1
    % At least 80% of global maximum
    peakIntervals1 = find(abs(crossCorr{1})>(0.8*max(abs(crossCorr{1}))));
    % We know that the training signal has more than 1000 samples
    peakCandidates1 = find((peakIntervals1(2:end)-peakIntervals1(1:end-1))>1000);
%     if ~isempty(peakCandidates1)
%         peakCandidates1 = peakCandidates1(1:end-1);  % Don't use the last one
%     end
    channelEstimate1 = zeros(nTxAntennas1, length(peakCandidates1));
    for i = 1:length(peakCandidates1)
        startIndex = peakIntervals1(peakCandidates1(i)+1);
        myStartIdx1(i) = startIndex;  %#ok<AGROW>
        % Focus on a small window that contains a peak. +/- 100 samples
        if startIndex > 100; windowRange = startIndex+(-100:100);
        else;                windowRange = 1:(startIndex+100);
        end

        dataWindow = [];
        for txID = 1:nTxAntennas1
            dataWindow = [dataWindow crossCorr{txID}(windowRange)];  %#ok<AGROW>
        end
        dataWindow = abs(dataWindow);

        % Find exact locations of the peaks (one in each of 4 cross correlations)
        [~,maxLoc] = max(dataWindow); % maxLoc is a 1x4 vector
        finalPeakLoc = startIndex-101;
            
        indexAll = finalPeakLoc + maxLoc + rxSigLen;
        indexMin = min(indexAll);
        indexMax = max(indexAll);

        if (indexMax-rxSigLen+L+400+payloadLen) <= rxSigLen % 400 samples between training signal and payload
            % Got a complete payload
            % Extract received training signal and payload
            receivedTrainingSig = rxSig((indexMin-rxSigLen+1):(indexMax-rxSigLen+L));
            payload = payload + ...
              rxSig((indexMax-rxSigLen+L+400+1):(indexMax-rxSigLen+L+400+payloadLen));
            numPayload = numPayload + 1;

            % Training signals from the 4 TX antennas may arrive at the RX at slightly different times.
            % Potentially off by 1 or 2 samples. Need to align Gold sequences before least squares fit.
            refSigShifted = zeros(L+indexMax-indexMin,nTxAntennas1);
            for txID = 1:nTxAntennas1
                refSigShifted((indexAll(txID)-indexMin)+(1:L),txID) = goldSeqRef(:,txID);
            end

            % Use least squares fit to estimate channel response
            channelEstimate1(:,i) = refSigShifted\receivedTrainingSig;
        else
            % Got an incomplete payload at the end of rxSig
            channelEstimate1(:,i) = [];
            break;
        end
    end
    channelEstimate1 = mean(channelEstimate1.');
else
    channelEstimate1 = [];
end

%% FOR TX2
% Find intervals containing peaks of cross correlation FOR CHANNEL 1 - TX2
% At least 80% of global maximum
if nTxAntennas2 > 0
    peakIntervals2 = find(abs(crossCorr{nTxAntennas1+1})>(0.8*max(abs(crossCorr{nTxAntennas1+1}))));
    % We know that the training signal has more than 1000 samples
    peakCandidates2 = find((peakIntervals2(2:end)-peakIntervals2(1:end-1))>1000);
%     if ~isempty(peakCandidates2)
%         peakCandidates2 = peakCandidates2(1:end-1);  % Don't use the last one
%     end
    channelEstimate2 = zeros(nTxAntennas2, length(peakCandidates2));
    for i = 1:length(peakCandidates2)
        startIndex = peakIntervals2(peakCandidates2(i)+1);
        if i == 1  && ~isempty(myStartIdx1)  && length(myStartIdx1)<5
            v = myStartIdx1 - startIndex;
            [~,idx] = min(abs(v));
            delaySamp = v(idx);
            if abs(delaySamp) < 2;    finalCorrectionTime = [0 0];
            elseif delaySamp > 33e3;  finalCorrectionTime = [0 0];
            elseif delaySamp > 0;     finalCorrectionTime = [0 delaySamp];
            else;                     finalCorrectionTime = [abs(delaySamp) 0];
            end
        end

        % Focus on a small window thatcd Desktop/BFInfocom/ contains a peak. +/- 100 samples
        if startIndex > 100; windowRange = startIndex+(-100:100);
        else;                windowRange = 1:(startIndex+100);
        end

        dataWindow = [];
        for txID = 1:nTxAntennas2
            dataWindow = [dataWindow crossCorr{nTxAntennas1+txID}(windowRange)];  %#ok<AGROW>
        end
        dataWindow = abs(dataWindow);

        % Find exact locations of the peaks (one in each of 4 cross correlations)
        [~,maxLoc] = max(dataWindow); % maxLoc is a 1x4 vector
        finalPeakLoc = startIndex-101;

        indexAll = finalPeakLoc + maxLoc + rxSigLen;
        indexMin = min(indexAll);
        indexMax = max(indexAll);

        if (indexMax-rxSigLen+L+400+payloadLen) <= rxSigLen % 400 samples between training signal and payload
            % Got a complete payload
            % Extract received training signal and payload
            receivedTrainingSig = rxSig((indexMin-rxSigLen+1):(indexMax-rxSigLen+L));
            payload = payload + ...
              rxSig((indexMax-rxSigLen+L+400+1):(indexMax-rxSigLen+L+400+payloadLen));
            numPayload = numPayload + 1;

            % Training signals from the 4 TX antennas may arrive at the RX at slightly different times.
            % Potentially off by 1 or 2 samples. Need to align Gold sequences before least squares fit.
            refSigShifted = zeros(L+indexMax-indexMin,nTxAntennas1);
            for txID = 1:nTxAntennas2
                refSigShifted((indexAll(txID)-indexMin)+(1:L),txID) = goldSeqRef(:,2+txID);
            end

            % Use least squares fit to estimate channel response
            channelEstimate2(:,i) = refSigShifted\receivedTrainingSig;
        else
            % Got an incomplete payload at the end of rxSig
            channelEstimate2(:,i) = [];
            break;
        end
    end
    channelEstimate2 = mean(channelEstimate2.');
else
    channelEstimate2 = [];
end

channelEstimate = [channelEstimate1 channelEstimate2];

% Evaluate time correction between transmitters
fprintf('System needs to be corrected by %d and %d samples\n',finalCorrectionTime);

if numPayload > 0
    % The transmitter sends the same payload continously.
    % Return the average of all detected payloads in rxSig.
    payload = payload/numPayload;
end
