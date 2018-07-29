clear all; close all; clc;

maxIter = 300;
tol = 20;
iters = (tol:1:maxIter);
antList = [1 2 4];
gainList = [20 10 0];
modList = [64 32 16 8 4 2];

%% GAIN
h =  findobj('type','figure');
n = length(h);
figure(n + 1); hold on;
nTxAntennasList = 4;
power = 0;
leg = cell(length(nTxAntennasList),1);
for nTxAntennas = nTxAntennasList
    dataFile = strcat('data/expOutdoor_',num2str(nTxAntennas),'Tx1Rx_TxHov_RxStat_',...
                      num2str(power),'dBm_1-5m.mat');
    load(dataFile);
    totGain = 0;
    for txID = 1:nTxAntennas
        t = abs(chTot(txID,1:maxIter));
        t(isnan(t))=0;  % replace nan values
        totGain = totGain + t;
    end
    plot((1:maxIter),totGain,'LineWidth',2);
    leg(nTxAntennas==nTxAntennasList) = strcat('NTx Antennas:',{' '},num2str(nTxAntennas));
end
figure(n + 1);
title('Total channel Gain','FontSize',14);
xlabel('Iteration','FontSize',14);
ylabel('Gain (linear)','FontSize',14);
hleg = legend(leg);
set(hleg,'FontSize',12);

h =  findobj('type','figure');  n = length(h);
nTxAntennas = 4;
power = 0;
dataFile = strcat('data/expOutdoor_',num2str(nTxAntennas),'Tx1Rx_TxStat_RxStat_',...
                      num2str(power),'dBm_1-5m.mat');
load(dataFile);
for idx = 1:nTxAntennas
    figure(n + 1);
    subplot(nTxAntennas,1,idx); hold on;
    covTot_imag = zeros(maxIter,1);
    covTot_real = zeros(maxIter,1);
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
    title(strcat('Ant id:',{' '},num2str(idx)),'FontSize',14);
%     ylim([0.999 1]);

    figure(n + 2);
    subplot(nTxAntennas,3,3*(idx-1) + 1); hold on; grid minor;
    % plot((1:maxIter),real(chRealTot),'Color','b');
    plot((1:maxIter),real(chTot(idx,1:maxIter)),'Color','r','LineWidth',2);
    title('Real','FontSize',12);
    xlabel('Iteration','FontSize',12);
    ylabel('Gain','FontSize',12);
    subplot(nTxAntennas,3,3*(idx-1) + 2); hold on; grid minor;
    % plot((1:maxIter),imag(chRealTot),'Color','b');
    plot((1:maxIter),imag(chTot(idx,1:maxIter)),'Color','r','LineWidth',2);
    title('Imaginary','FontSize',12);
    xlabel('Iteration','FontSize',12);
    ylabel('Gain','FontSize',12);
    subplot(nTxAntennas,3,3*(idx-1) + 3); hold on; grid minor;
    plot((1:maxIter),abs(chTot(idx,1:maxIter)),'Color','r','LineWidth',2);
    title('Absolute Gain','FontSize',12);
    xlabel('Iteration','FontSize',12);
    ylabel('Gain','FontSize',12);
end

%% ECDF - over antenna selection
h =  findobj('type','figure');  n = length(h);
figure(n + 1); hold on;
power = 20;
antList = [1 2 4];
myModList = [64];
powIdx = find(gainList==power);
leg = cell(length(myModList),length(antList));
legIdx = 1;
for ant = antList
    antIdx = find(ant==antList);
    for mod = myModList
        modIdx = find(mod==myModList);
        [Y,X] = ecdf(BER(iters,modIdx,antIdx,powIdx));
        plot(X.*100,Y,'LineWidth',2);
        leg(legIdx) = strcat('Ntx:',{' '},num2str(ant),{' '},...
                                     'Modulation:',{' '},num2str(mod),'QAM');
        legIdx = legIdx + 1;
    end
end
hleg = legend(leg);
set(hleg,'FontSize',12);
xlabel('BER(%)','FontSize',14);
ylabel('F(x)','FontSize',14);
title('BER distribution','FontSize',14);

%%
power = 0;
nTxAntennasList = [1 2 4];
cap_max = zeros(length(nTxAntennasList),1);
cap_achieved = zeros(length(nTxAntennasList),1);
for txID = 1:length(nTxAntennasList)
    nTxAntennas = nTxAntennasList(txID);
    dataFile = strcat('data/expOutdoor_',num2str(nTxAntennas),'Tx1Rx_TxStat_RxStat_',...
                          num2str(power),'dBm_1-5m.mat');
    load(dataFile,'chTot');
    power_lin = db2pow(power-10);
    channel = normalize(chTot(:,tol:maxIter-1));  % Take one sample at the end
    channel_est = normalize(chTot(:,tol+1:maxIter));  % Shift channel by 1 iterations
    ch1 = 0;  ch2 = 0;
    for k = 1:nTxAntennasList(txID)
        ch1 = ch1 + abs(channel(k,:) .* conj(channel(k,:)));
        ch2 = ch2 + abs(channel_est(k,:) .* conj(channel(k,:)));
    end
    cap_max(txID) = log10(1 + power_lin.*ch1);
    cap_achieved(txID) = log10(1 + power_lin.*ch2);
end
h =  findobj('type','figure');  n = length(h);
figure(n + 1); hold on;
plot(nTxAntennasList,cap_max,'LineWidth',2);
plot(nTxAntennasList,cap_achieved,'LineWidth',2);

%%
% load('BERTot_file','BER');
load('BERTot_Rx_file','BER_RxHov');
BER = BER_RxHov;

%% Comparison amongst antenna configurations (mean)

figure(10); hold on;
myModList = [64 32 16 8 4 2];

for mod = myModList
    x_1Tx = BER(tol:maxIter,mod==modList,1,1);
    x_2Tx = BER(tol:maxIter,mod==modList,2,1);
    x_4Tx = BER(tol:maxIter,mod==modList,3,1);

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

    data_tot = [mean(x_1Tx) mean(x_2Tx) mean(x_4Tx)].*100;
    error_tot = [error_1Tx error_2Tx error_4Tx].*100;
    errorbar(antList,data_tot,error_tot,'-s','LineWidth',2);
    leg(mod==myModList) = strcat('Modulation:',{' '},num2str(mod));
end
hleg = legend(leg);
set(hleg,'FontSize',12);
xlabel('Number of antennas','FontSize',14);
ylabel('BER(%)','FontSize',14);
title('BER improvement with Beamforming','FontSize',14);