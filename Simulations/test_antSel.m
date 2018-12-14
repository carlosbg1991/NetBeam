clear all; clear classes; close all; clc;

N = 20;
MList = (1:1:10);
Niter = 10;
% myAlphas = [0.1 0.3 0.5 0.9 1].';
myAlphas = 1;

% Generate channel
chMaxTot = (random('Rician',0,1,max(MList),N)+1i*random('Rician',0,1,max(MList),N));

%% TEST 1
% % Determine node replicabiliy
% maxAssignationList = 1:1:N;  % Cover all the possibilities
% Csum = zeros(length(maxAssignationList),length(MList));
% lenAssig = zeros(max(MList),length(MList));
% for idxM = 1:length(MList)
%     for iter = 1:Niter
%          
%         M = MList(idxM);
%         
%         % Generate channel
%         chMax = chMaxTot(1:M,:);
%         
%         % Generate path loss
%         pos = randi(length(myAlphas),M,1);
%         alphas = myAlphas(pos);
%         chMax = alphas.*chMax;
% 
%         % Cost Matrix as the Input for the Hungarian Algorithm
%         chMaxGain = abs(chMax);
%         myMax = max(chMaxGain,[],'all');
%         chMaxGain_norm = chMaxGain./myMax + 1e-5;
%         W = 1 - chMaxGain_norm;  % Initial channel gain (the higher, the better)
% 
%         for maxAssignation = maxAssignationList
% 
%             % Replicate transmitter nodes
%             W_Mod = []; pVec = [];
%             % W_Mod = zeros(replica,N);
%             for txID = 1:M
%                 W1 = repmat(W(txID,:),maxAssignation,1);
%                 W_Mod = [W_Mod ; W1];                                 %#ok<AGROW>
%             end
% 
%             % Hungarian Algorithm
%             [ Assign , ~ ] = munkres(W_Mod);
% 
%             % Upper-bound Capacity 
%             C = zeros(1,M);
%             for txID = 1:M
%                 % Indices
%                 myIni = (txID-1)*(maxAssignation) + 1;
%                 myEnd = txID*maxAssignation;
% 
%                 % transmitters assigned to user txID
%                 myAssign = Assign(myIni:myEnd);
%                 myAssign(myAssign==0) = [];
%                 lenAssig(txID,idxM) = length(myAssign);
% 
%                 myChannel = chMaxGain(txID,myAssign);
%                 C(txID) = sum(myChannel);
%             end
% 
%             Csum(maxAssignation==maxAssignationList,idxM) = sum(C);
%             fprintf('%d - %d ---> %.4f\n',N,M,Csum(maxAssignation==maxAssignationList));
%             fprintf('%d | ',lenAssig);
%             fprintf('\n')
% 
%         end
%         
%         figure(1); hold on;
%         plot(maxAssignationList,Csum(:,idxM));
%         maxCumSum(idxM) = max(Csum(:,idxM));
%     end
% end
% 
% figure(2);
% plot(MList,maxCumSum);


%% TEST 2: Just evaluate it using the maximum: N
maxAssignation = N;
Csum = zeros(length(MList),Niter);
lenAssig = zeros(max(MList),length(MList));

for iter = 1:Niter

    % Generate channel
    chMaxTot = (random('Rician',0,1,max(MList),N)+1i*random('Rician',0,1,max(MList),N));

    for idxM = 1:length(MList)

            M = MList(idxM);

            % Generate channel
            chMax = chMaxTot(1:M,:);

            % Generate path loss
            pos = randi(length(myAlphas),M,1);
            alphas = myAlphas(pos);
            chMax = alphas.*chMax;

            % Cost Matrix as the Input for the Hungarian Algorithm
            chMaxGain = abs(chMax);
            myMax = max(chMaxGain,[],'all');
            chMaxGain_norm = chMaxGain./myMax + 1e-5;
            W = 1 - chMaxGain_norm;  % Initial channel gain (the higher, the better)

            % Replicate transmitter nodes
            W_Mod = []; pVec = [];
            % W_Mod = zeros(replica,N);
            for txID = 1:M
                W1 = repmat(W(txID,:),maxAssignation,1);
                W_Mod = [W_Mod ; W1];                                 %#ok<AGROW>
            end

            % Hungarian Algorithm
            [ Assign , ~ ] = munkres(W_Mod);

            % Upper-bound Capacity 
            C = zeros(1,M);
            for txID = 1:M
                % Indices
                myIni = (txID-1)*(maxAssignation) + 1;
                myEnd = txID*maxAssignation;

                % transmitters assigned to user txID
                myAssign = Assign(myIni:myEnd);
                myAssign(myAssign==0) = [];
                lenAssig(txID,idxM) = length(myAssign);

                myChannel = chMaxGain(txID,myAssign);
                C(txID) = sum(myChannel);
            end

            Csum(idxM,iter) = sum(C);
            fprintf('%d - %d ---> %.4f\n',N,M,Csum(idxM));
            fprintf('%d | ',lenAssig);
            fprintf('\n')
    end

end

Csum_av = mean(Csum,2);

figure(1); hold on;
plot(MList,Csum);
figure(2); hold on;
plot(MList,Csum_av);

