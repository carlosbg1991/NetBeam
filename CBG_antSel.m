function [TxAssignedToRx] = CBG_antSel(chMax,K)

% chMax: Matrix of size M by N, where N is the number of transmitters and M
% the number of receivers.

% Cost Matrix as the Input for the Hungarian Algorithm
W = chMax;  % Initial channel gain (the higher, the better)
N = size(chMax,2);  % Number of transmitters
M = size(chMax,1);  % Number of receivers

% Initialization of the Modified Hungarian Algorithm
P = ceil(ones(size(W,1),size(W,2)) - (ones(size(W,1),size(W,2)) - W).^50);
t = ones(1,size(P,1))*P;
dummy_relays = t;  % Check Sec. IV for proof

% Dummy Relays insertion
W_Mod = []; pVec = [];
for col = 1:size(W,2)-1    %The last one belongs to the AP
    W1 = repmat(W(:,col),1, min(dummy_relays(col),K));
    W_Mod = [W_Mod W1];                                   %#ok<AGROW>
    pV1 = repmat(col, 1, min(dummy_relays(col),K));
    pVec = [pVec pV1];                                    %#ok<AGROW>
end

% Hungarian Algorithm
W_Mod_reverse = ones(size(W_Mod,1),size(W_Mod,2)) - W_Mod;
[ Assign , ~ ] = munkres(W_Mod_reverse);
TxAssignedToRx = pVec(Assign);

% EOF