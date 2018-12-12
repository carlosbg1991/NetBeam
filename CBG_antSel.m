function [finalAssign,SNR] = CBG_antSel(chMax,K)
% CBG_antSel - returns the antenna allocation given the channel matrix and
% the maximum channel replicas, to balance out between maximum sum capacity
% and equity across users.
%
% Syntax:  [Assign,SNR] = CBG_antSel(chMax,K)
%
% Inputs:
%    chMax - Matrix of size M by N, where M the number of receivers.N is
%            the number of transmitters and
%    K     - Maximum number of transmitters replicas
%
% Outputs:
%    finalAssign - Cell, distribution of transmitters across users. For
%             instance: Assign{2} = [28 37 44 12] shows the IDs of the
%             transmitters assigned to user 2.
%    SNR -    Scalar, contains the maximum SNR achieved per user
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

% Parse inputs
N = size(chMax,2);  % Number of transmitters
M = size(chMax,1);  % Number of receivers

% Cost Matrix as the Input for the Hungarian Algorithm

% Cost Matrix as the Input for the Hungarian Algorithm
chMaxGain = abs(chMax);
myMax = max(chMaxGain,[],'all');
chMaxGain_norm = chMaxGain./myMax;
W = 1 - chMaxGain_norm;  % Initial channel gain (the higher, the better)

% Replicate transmitter nodes
W_Mod = [];
for txID = 1:M
    W1 = repmat(W(txID,:),K,1);
    W_Mod = [W_Mod ; W1];   %#ok<AGROW>
end

% Hungarian Algorithm
[ Assign , ~ ] = munkres(W_Mod);

% Upper-bound Capacity 
SNR = zeros(1,M);
finalAssign = cell(1,M);
for txID = 1:M
    % Indices
    myIni = (txID-1)*K + 1;
    myEnd = txID*K;

    % transmitters assigned to user txID
    myAssign = Assign(myIni:myEnd);
    myAssign(myAssign==0) = [];
    finalAssign{txID} = myAssign;

    myChannel = chMaxGain(txID,myAssign);
    SNR(txID) = sum(myChannel);
end

% Csum = sum(SNR);
% fprintf('%d TX  -  %d RX  -> max capacity %.4f\n',N,M,Csum);


% EOF