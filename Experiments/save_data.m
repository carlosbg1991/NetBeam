%%
maxIter = 300;
tol = 20;
iters = (tol:1:maxIter);
antList = [1 2 4];
gainList = [20 10 0];
modList = [64 32 16 8 4 2];

%% STATIC and OUTDOOR
% BER = zeros(300,length(modList),length(antList),length(gainList));
%
% data = load('data/expOutdoor_1Tx1Rx_TxStat_RxStat_-20dBm_1-5m.mat');
% BER(:,:,1,1) = data.BER(1:maxIter,:);
% data_1Tx_neg10 = load('data/expOutdoor_1Tx1Rx_TxStat_RxStat_-10dBm_1-5m.mat');
% BER(:,:,1,2) = data.BER(1:maxIter,:);
% data_1Tx_0 = load('data/expOutdoor_1Tx1Rx_TxStat_RxStat_0dBm_1-5m.mat');
% BER(:,:,1,3) = data.BER(1:maxIter,:);
% % 
% data = load('data/expOutdoor_2Tx1Rx_TxStat_RxStat_-20dBm_1-5m.mat');
% BER(:,:,2,1) = data.BER(1:maxIter,:);
% data_2Tx_neg10 = load('data/expOutdoor_2Tx1Rx_TxStat_RxStat_-10dBm_1-5m.mat');
% BER(:,:,2,2) =data.BER(1:maxIter,:);
% data_2Tx_0 = load('data/expOutdoor_2Tx1Rx_TxStat_RxStat_0dBm_1-5m.mat');
% BER(:,:,2,2) = data.BER(1:maxIter,:);
% % 
% data = load('data/expOutdoor_4Tx1Rx_TxStat_RxStat_-20dBm_1-5m.mat');
% BER(:,:,3,1) = data.BER(1:maxIter,:);
% data_4Tx_neg10 = load('data/expOutdoor_4Tx1Rx_TxStat_RxStat_-10dBm_1-5m.mat');
% BER(:,:,3,2) = data.BER(1:maxIter,:);
% data_4Tx_0 = load('data/expOutdoor_4Tx1Rx_TxStat_RxStat_0dBm_1-5m.mat');
% BER(:,:,3,3) = data.BER(1:maxIter,:);
% save('BERTot_file','BER');

%% RX-HOV and OUTDOOR
% BER_RxHov = zeros(300,length(modList),length(antList),length(gainList));
% 
% data = load('data/expOutdoor_1Tx1Rx_TxStat_RxHov_-20dBm_1-5m.mat');
% BER_RxHov(:,:,1,1) = data.BER(1:maxIter,:);
% data_1Tx_neg10 = load('data/expOutdoor_1Tx1Rx_TxStat_RxHov_-10dBm_1-5m.mat');
% BER_RxHov(:,:,1,2) = data.BER(1:maxIter,:);
% data_1Tx_0 = load('data/expOutdoor_1Tx1Rx_TxStat_RxHov_0dBm_1-5m.mat');
% BER_RxHov(:,:,1,3) = data.BER(1:maxIter,:);
% % 
% data = load('data/expOutdoor_2Tx1Rx_TxStat_RxHov_-20dBm_1-5m.mat');
% BER_RxHov(:,:,2,1) = data.BER(1:maxIter,:);
% data_2Tx_neg10 = load('data/expOutdoor_2Tx1Rx_TxStat_RxHov_-10dBm_1-5m.mat');
% BER_RxHov(:,:,2,2) =data.BER(1:maxIter,:);
% data_2Tx_0 = load('data/expOutdoor_2Tx1Rx_TxStat_RxHov_0dBm_1-5m.mat');
% BER_RxHov(:,:,2,2) = data.BER(1:maxIter,:);
% % 
% data = load('data/expOutdoor_4Tx1Rx_TxStat_RxHov_-20dBm_1-5m.mat');
% BER_RxHov(:,:,3,1) = data.BER(1:maxIter,:);
% data_4Tx_neg10 = load('data/expOutdoor_4Tx1Rx_TxStat_RxHov_-10dBm_1-5m.mat');
% BER_RxHov(:,:,3,2) = data.BER(1:maxIter,:);
% data_4Tx_0 = load('data/expOutdoor_4Tx1Rx_TxStat_RxHov_0dBm_1-5m.mat');
% BER_RxHov(:,:,3,3) = data.BER(1:maxIter,:);
% save('BERTot_Rx_file','BER_RxHov');

%% TX-HOV and OUTDOOR
% BER_TxHov = zeros(300,length(modList),length(antList),length(gainList));
% 
% data = load('data/expOutdoor_1Tx1Rx_TxHov_RxStat_-20dBm_1-5m.mat');
% BER_TxHov(:,:,1,1) = data.BER(1:maxIter,:);
% data_1Tx_neg10 = load('data/expOutdoor_1Tx1Rx_TxHov_RxStat_-10dBm_1-5m.mat');
% BER_TxHov(:,:,1,2) = data.BER(1:maxIter,:);
% data_1Tx_0 = load('data/expOutdoor_1Tx1Rx_TxHov_RxStat_0dBm_1-5m.mat');
% BER_TxHov(:,:,1,3) = data.BER(1:maxIter,:);
% % 
% data = load('data/expOutdoor_2Tx1Rx_TxHov_RxStat_-20dBm_1-5m.mat');
% BER_TxHov(:,:,2,1) = data.BER(1:maxIter,:);
% data_2Tx_neg10 = load('data/expOutdoor_2Tx1Rx_TxHov_RxStat_-10dBm_1-5m.mat');
% BER_TxHov(:,:,2,2) =data.BER(1:maxIter,:);
% data_2Tx_0 = load('data/expOutdoor_2Tx1Rx_TxHov_RxStat_0dBm_1-5m.mat');
% BER_TxHov(:,:,2,2) = data.BER(1:maxIter,:);
% % 
% data = load('data/expOutdoor_4Tx1Rx_TxHov_RxStat_-20dBm_1-5m.mat');
% BER_TxHov(:,:,3,1) = data.BER(1:maxIter,:);
% data_4Tx_neg10 = load('data/expOutdoor_4Tx1Rx_TxHov_RxStat_-10dBm_1-5m.mat');
% BER_TxHov(:,:,3,2) = data.BER(1:maxIter,:);
% data_4Tx_0 = load('data/expOutdoor_4Tx1Rx_TxHov_RxStat_0dBm_1-5m.mat');
% BER_TxHov(:,:,3,3) = data.BER(1:maxIter,:);
% save('BERTot_Tx_file','BER_TxHov');

%% TX-HOV and KRI
% The only gain used throughout the 4 configurations is gain=10dB
gainList = 15;
BER_TxHov_KRI = zeros(300,length(modList),length(antList),length(gainList));

data = load('data/expKRI_1Tx1Rx_TxHov_RxStat_15dB_5m.mat');
BER_TxHov_KRI(:,:,1,1) = data.BER(1:maxIter,:);
data = load('data/expKRI_2Tx1Rx_TxHov_RxStat_15dB_5m.mat');
BER_TxHov_KRI(:,:,2,1) = data.BER(1:maxIter,:);
data = load('data/expKRI_4Tx1Rx_TxHov_RxStat_15dB_5m.mat');
BER_TxHov_KRI(:,:,2,1) = data.BER(1:maxIter,:);
save('BERTot_Tx_KRI_file','BER_TxHov_KRI');