maxIter = 300;
tol = 20;
iters = (tol:1:maxIter);
antList = [1 2 4];
gainList = [20 10 0];
modList = [64 32 16 8 4 2];
BER = zeros(300,length(modList),length(antList),length(gainList));

data = load('data/expOutdoor_1Tx1Rx_TxStat_RxStat_-20dBm_1-5m.mat');
BER(:,:,1,1) = data.BER(1:maxIter,:);
data_1Tx_neg10 = load('data/expOutdoor_1Tx1Rx_TxStat_RxStat_-10dBm_1-5m.mat');
BER(:,:,1,2) = data.BER(1:maxIter,:);
data_1Tx_0 = load('data/expOutdoor_1Tx1Rx_TxStat_RxStat_0dBm_1-5m.mat');
BER(:,:,1,3) = data.BER(1:maxIter,:);
% 
data = load('data/expOutdoor_2Tx1Rx_TxStat_RxStat_-20dBm_1-5m.mat');
BER(:,:,2,1) = data.BER(1:maxIter,:);
data_2Tx_neg10 = load('data/expOutdoor_2Tx1Rx_TxStat_RxStat_-10dBm_1-5m.mat');
BER(:,:,2,2) =data.BER(1:maxIter,:);
data_2Tx_0 = load('data/expOutdoor_2Tx1Rx_TxStat_RxStat_0dBm_1-5m.mat');
BER(:,:,2,2) = data.BER(1:maxIter,:);
% 
data = load('data/expOutdoor_4Tx1Rx_TxStat_RxStat_-20dBm_1-5m.mat');
BER(:,:,3,1) = data.BER(1:maxIter,:);
data_4Tx_neg10 = load('data/expOutdoor_4Tx1Rx_TxStat_RxStat_-10dBm_1-5m.mat');
BER(:,:,3,2) = data.BER(1:maxIter,:);
data_4Tx_0 = load('data/expOutdoor_4Tx1Rx_TxStat_RxStat_0dBm_1-5m.mat');
BER(:,:,3,3) = data.BER(1:maxIter,:);

%% ECDF - over antenna selection
figure(1); hold on;
power = 10;
antList = [1 2 4];
myModList = [32];

powIdx = find(gainList==power);
for ant = antList
    antIdx = find(ant==antList);
    for mod = myModList
        modIdx = find(mod==myModList);
        ecdf(BER(iters,modIdx,antIdx,powIdx));
    end
end
legend('64QAM','32QAM','64QAM','32QAM');

%% ECDF - over antenna selection
figure(1); hold on;
ant = 1;
power = 0;
myModList = [64 32 16];

antIdx = find(antList==ant);
powIdx = find(gainList==power);
for mod = myModList
    modIdx = find(mod==myModList);
    ecdf(BER(iters,modIdx,antIdx,powIdx));
end
legend('64QAM','32QAM','16QAM');

%% Comparison amongst antenna configurations (mean)
lastIter = 300;
tol = 20;

x_1Tx = data_1Tx_20.BER(tol:lastIter,1);
x_2Tx = data_2Tx_20.BER(tol:lastIter,1);
x_4Tx = data_4Tx_20.BER(tol:lastIter,1);

% Error 1 UAV
SEM = std(x_1Tx)/sqrt(length(x_1Tx));           % Standard Error
ts = tinv([0.025  0.975],length(x_1Tx)-1);      % T-Score
CI = mean(x_1Tx) + ts*SEM;                      % Confidence Intervals
error_1Tx = abs(CI(1) - mean(x_1Tx));
% Error 2 UAV's
SEM = std(x_2Tx)/sqrt(length(x_2Tx));               % Standard Error
ts = tinv([0.025  0.975],length(x_2Tx)-1);      % T-Score
CI = mean(x_2Tx) + ts*SEM;                      % Confidence Intervals
error_2Tx = abs(CI(1) - mean(x_2Tx));
% Error 4 UAV's
SEM = std(x_4Tx)/sqrt(length(x_4Tx));               % Standard Error
ts = tinv([0.025  0.975],length(x_4Tx)-1);      % T-Score
CI = mean(x_4Tx) + ts*SEM;                      % Confidence Intervals
error_4Tx = abs(CI(1) - mean(x_4Tx));


figure(10); hold on;
data_tot = [mean(x_1Tx) mean(x_2Tx) mean(x_4Tx)];
error_tot = [error_1Tx error_2Tx error_4Tx];
errorbar((1:length(data_tot)),data_tot,error_tot,'-s');