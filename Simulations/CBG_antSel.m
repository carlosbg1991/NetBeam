function [finalAssign,SNR,assignation] = CBG_antSel(chMax,SNRdemand,policy)
% CBG_antSel - returns the antenna allocation given the channel matrix and
% the maximum channel replicas, to balance out between maximum sum capacity
% and equity across users.
%
% Syntax:  [Assign,SNR] = CBG_antSel(chMax,K)
%
% Inputs:
%    chMax - Matrix of size M by N, where M the number of receivers.N is
%            the number of transmitters and M the number of receivers.
%    SNRdemand - Vector, SNR demanded
%    policy - Selection policy: 'optimum', 'random', 'greedy',
%             'greedy-equtv'
%
% Outputs:
%    finalAssign - Cell, distribution of transmitters across users. For
%             instance: Assign{2} = [28 37 44 12] shows the IDs of the
%             transmitters assigned to user 2.
%    SNR -    Scalar, contains the maximum SNR achieved per user
%    channel - Matrix of size M by N, where M the number of receivers.N is
%              the number of transmitters and M the number of receivers.
%
% Example: 
%    To-do
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: To-do

%------------- BEGIN CODE --------------

% Priority list
SNRdemand = SNRdemand.';  % Need to correct
[~,priority] = sort(SNRdemand,'descend');

% Parse inputs
N = size(chMax,2);  % Number of transmitters
M = size(chMax,1);  % Number of receivers

% Decide how many antennas we are going to assign
SNRdemandNorm = SNRdemand./sum(SNRdemand);
K = floor(SNRdemandNorm*N);
% Add unasigned transmitter to user with highest demands
if sum(K) < N
    K(priority(1)) = K(priority(1)) + 1;
end

% Cost Matrix as the Input for the Hungarian Algorithm
chMaxGain = abs(chMax);
myMax = max(chMaxGain,[],'all');
chMaxGain_norm = chMaxGain./myMax;
W = 1 - chMaxGain_norm;  % Initial channel gain (the higher, the better)

finalAssign = cell(1,M);
SNR = zeros(1,M);
if strcmp(policy,'optimum')   % Optimum transmitter allocation following the modified Hungarian
    % Replicate transmitter nodes
    W_Mod = [];
    for rxID = 1:M
        W1 = repmat(W(rxID,:),K(rxID),1);
        W_Mod = [W_Mod ; W1];   %#ok<AGROW>
    end

    % Hungarian Algorithm
    [ Assign , ~ ] = munkres(W_Mod);

    myEnd = 0;
    for rxID = 1:M
        % Indices
        myIni = myEnd + 1;
        myEnd = myIni + K(rxID) - 1;

        % transmitters assigned to user txID
        myAssign = Assign(myIni:myEnd);
        myAssign(myAssign==0) = [];
        finalAssign{rxID} = myAssign;

        myChannel = chMaxGain(rxID,myAssign);
        SNR(rxID) = sum(myChannel);
    end
elseif strcmp(policy,'random')   % Random policy selection
    totPoss = 1:1:N;
    for rxID = priority
        [~,idxShuffle] = datasample(totPoss,...
                            min(length(totPoss),ceil(N/M)),'Replace',false);
        
        finalAssign{rxID} = totPoss(idxShuffle);
        totPoss(idxShuffle) = [];
    end
elseif strcmp(policy,'greedy')
    totPoss = 1:1:N;  % List of Tx to be assigned to Rx
    while(~isempty(totPoss))
        for rxID = priority
            [~,greedyIdx] = sort(chMaxGain(rxID,totPoss),'descend');
            mySelect = greedyIdx(1:ceil(N/M));
            finalAssign{rxID} = [finalAssign{rxID} totPoss(mySelect)];
            totPoss(mySelect) = [];
        end
    end
elseif strcmp(policy,'greedy-equtv')
    totPoss = 1:1:N;  % List of Tx to be assigned to Rx
    while(~isempty(totPoss))
        for rxID = priority
            [~,greedyIdx] = max(chMaxGain(rxID,totPoss));
            finalAssign{rxID} = [finalAssign{rxID} totPoss(greedyIdx)];
            totPoss(greedyIdx) = [];
        end
    end
end

% Compute achieved SNR
for rxID = 1:M
    myChannel = chMaxGain(rxID,finalAssign{rxID});
    SNR(rxID) = sum(myChannel);
end

%% Compute modified channel matrix
assignation = zeros(N,M);
for rxID = 1:M
	assignation(finalAssign{rxID},rxID) = 1;
end

% Generation of samples using Raylegh and Rician
% N = 1e6; close all
% chMax = (raylrnd(1,[N 1])+1i*raylrnd(1,[N 1]));
% histogram(abs(chMax),100); hold on;
% chMax = (random('Rician',0,2,N,1)+1i*random('Rician',0,2,N,1));
% histogram(abs(chMax),100)
% grid minor
% EOF