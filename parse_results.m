clear all; close all; clc;

maxIter = 300;
tol = 20;
iters = (tol:1:maxIter);
antList = [1 2 4];
gainList = [20 10 0];
modList = [64 32 16 8 4 2];

%% BER
% h =  findobj('type','figure');
% n = length(h);
% figure(n + 1); 
nTxAnt = 4;
power = 15;
dataFile = strcat('data/expKRI_',num2str(nTxAnt),'Tx1Rx_TxHov_RxStat_',...
                      num2str(power),'dB_5m.mat');
% dataFile = strcat('data/expKRI_',num2str(nTxAntennas),'Tx1Rx_TxStat_RxStat_',...
%                       num2str(power),'dB_5m_Propsrunning.mat');
load(dataFile);
maxIter = 200;
subplot(611); hold on; plot(abs(BER(1:maxIter,1)),'LineWidth',1.3); legend('64-QAM'); ylabel('BER'); ylim([0 1]);
subplot(612); hold on; plot(abs(BER(1:maxIter,2)),'LineWidth',1.3); legend('32-QAM'); ylabel('BER'); ylim([0 1]);
subplot(613); hold on; plot(abs(BER(1:maxIter,3)),'LineWidth',1.3); legend('16-QAM'); ylabel('BER'); ylim([0 1]);
subplot(614); hold on; plot(abs(BER(1:maxIter,4)),'LineWidth',1.3); legend('8-QAM'); ylabel('BER'); ylim([0 1]);
subplot(615); hold on; plot(abs(BER(1:maxIter,5)),'LineWidth',1.3); legend('QPSK'); ylabel('BER'); ylim([0 1]);
subplot(616); hold on; plot(abs(BER(1:maxIter,6)),'LineWidth',1.3); legend('BPSK'); ylabel('BER'); ylim([0 1]);

%% GAIN
h =  findobj('type','figure');
n = length(h);
figure(n + 1); hold on;
nTxAntennasList = [1 2 4];
power = 15;
leg = cell(length(nTxAntennasList),1);
for nTxAnt = nTxAntennasList
%     dataFile = strcat('data/expOutdoor_',num2str(nTxAntennas),'Tx1Rx_TxHov_RxStat_',...
%                       num2str(power),'dBm_1-5m.mat');
    dataFile = strcat('data/expKRI_',num2str(nTxAnt),'Tx1Rx_TxHov_RxStat_',...
                      num2str(power),'dB_5m.mat');
    load(dataFile);
    maxIter = 300;
    totGain = 0;
    for txID = 1:nTxAnt
        t = abs(chTot(txID,1:maxIter));
        t(isnan(t))=0;  % replace nan values
        totGain = totGain + t;
    end
    plot((1:maxIter),totGain,'LineWidth',2);
    leg(nTxAnt==nTxAntennasList) = strcat('NTx Antennas:',{' '},num2str(nTxAnt));
end
figure(n + 1);
title('Total channel Gain - UAVs flying','FontSize',14);
xlabel('Iteration','FontSize',14);
ylabel('Gain (linear)','FontSize',14);
hleg = legend(leg);
set(hleg,'FontSize',12);


%% CHANNEL
h =  findobj('type','figure');  n = length(h);
nTxAnt = 4;
power = 15;
% dataFile = strcat('data/expOutdoor_',num2str(nTxAntennas),'Tx1Rx_TxStat_RxStat_',...
%                       num2str(power),'dBm_1-5m.mat');
dataFile = strcat('data/expKRI_',num2str(nTxAnt),'Tx1Rx_TxHov_RxStat_',...
                      num2str(power),'dB_5m.mat');
load(dataFile);
chTot_Hov = chTot;
dataFile = strcat('data/expKRI_',num2str(nTxAnt),'Tx1Rx_TxStat_RxStat_',...
                      num2str(power),'dB_5m_Propsrunning.mat');
load(dataFile);
chTot_Stat = chTot;
maxIter = 300;
for idx = 1:nTxAnt
    figure(n + 1);
    subplot(nTxAnt,1,idx); hold on;
    covTot_imag = zeros(maxIter,1);
    covTot_real = zeros(maxIter,1);
    re = real(chTot(idx,:));
    im = imag(chTot(idx,:));
    re(isnan(re))=0;
    window_size = 1;
    for k = 1:1:maxIter-window_size
        covTot_imag(k) = cov(im(k:k+window_size));
        covTot_real(k) = cov(re(k:k+window_size));
    end
    plot(abs(1-covTot_imag),'LineWidth',2);
    plot(abs(1-covTot_real),'LineWidth',2);
    legend('Imaginary','Real');
    title(strcat('Ant id:',{' '},num2str(idx)),'FontSize',14);
%     ylim([0.999 1]);

    figure(n + 2);
    subplot(nTxAnt,3,3*(idx-1) + 1); hold on; grid minor;
    plot((1:maxIter),real(chTot_Hov(idx,1:maxIter)),'Color','r','LineWidth',2);
    plot((1:maxIter),real(chTot_Stat(idx,1:maxIter)),'Color','b','LineWidth',2);
    title('Real','FontSize',12);
    xlabel('Iteration','FontSize',12);
    ylabel('Gain','FontSize',12);
    subplot(nTxAnt,3,3*(idx-1) + 2); hold on; grid minor;
    plot((1:maxIter),imag(chTot_Hov(idx,1:maxIter)),'Color','r','LineWidth',2);
    plot((1:maxIter),imag(chTot_Stat(idx,1:maxIter)),'Color','b','LineWidth',2);
    title('Imaginary','FontSize',12);
    xlabel('Iteration','FontSize',12);
    ylabel('Gain','FontSize',12);
    subplot(nTxAnt,3,3*(idx-1) + 3); hold on; grid minor;
    plot((1:maxIter),abs(chTot_Hov(idx,1:maxIter)),'Color','r','LineWidth',2);
    plot((1:maxIter),abs(chTot_Stat(idx,1:maxIter)),'Color','b','LineWidth',2);
    title('Absolute Gain','FontSize',12);
    xlabel('Iteration','FontSize',12);
    ylabel('Gain','FontSize',12);
end
figure(n+2);
leg = legend('Hovering','Static');
set(leg,'FontSize',12);

%%
load('BERTot_file','BER');
% load('BERTot_Rx_file','BER_RxHov');
% BER = BER_RxHov;
% load('BERTot_Tx_KRI_file','BER_TxHov_KRI');
% BER = BER_TxHov_KRI;

%% ECDF - over antenna selection
gainList = 15;
power = 15;
colorList = {[0 0.447 0.741], [0.85 0.325 0.098], [0.929 0.694 0.125]};
LineStyleList = {'-','-.','--',':'};
h =  findobj('type','figure');  n = length(h);
figure(n + 1); hold on;
antList = [1 2 4];
myModList = [64 32 16];
powIdx = find(gainList==power);
leg = cell(length(myModList)*length(antList),1);
legIdx = 1;
for ant = antList
    antIdx = find(ant==antList);
    for mod = myModList
        modIdx = find(mod==myModList);
        [Y,X] = ecdf(BER(tol:maxIter,modIdx,antIdx,powIdx));
        plot(X.*100,Y,'LineWidth',2,'Color',colorList{modIdx},'LineStyle',LineStyleList{antIdx});
        leg(legIdx) = strcat('Ntx:',{' '},num2str(ant),{' '},...
                                     'Modulation:',{' '},num2str(mod),'QAM');
        legIdx = legIdx + 1;
    end
end
hleg = legend(leg);
set(hleg,'FontSize',12);
xlabel('BER(%)','FontSize',14);
ylabel('F(x)','FontSize',14);
title('BER distribution - UAVs flying','FontSize',14);

%%
load('BERTot_file','BER');
% load('BERTot_Rx_file','BER_RxHov');
% BER = BER_RxHov;
% load('BERTot_Tx_KRI_file','BER_TxHov_KRI');
% BER = BER_TxHov_KRI;

%% CAPACITY (Theoretical)
% h =  findobj('type','figure');  n = length(h);
% figure(n + 1); hold on;
% power = 0;
power = 15;
nTxAntennasList = [1 2 4];
% ha = tight_subplot(length(nTxAntennasList),1,[.01 .01],[.05 .01],[.1 .3]);
ha = tight_subplot(length(nTxAntennasList),1,[.05 .05],[.05 .01],[.05 .18]);
cap_max = zeros(maxIter-tol,length(nTxAntennasList));
cap_achieved = zeros(maxIter-tol,length(nTxAntennasList));
for txID = 1:length(nTxAntennasList)
    nTxAnt = nTxAntennasList(txID);
%     dataFile = strcat('data/expOutdoor_',num2str(nTxAnt),'Tx1Rx_TxStat_RxStat_',...
%                           num2str(power),'dBm_1-5m.mat');
    dataFile = strcat('data/ExpKRI_',num2str(nTxAnt),'Tx1Rx_TxHov_RxStat_',...
                          num2str(power),'dB_5m.mat');
    load(dataFile,'chTot');
    power_lin = db2pow(power-10);
    ch1 = 0;  ch2 = 0;
    for k = 1:nTxAntennasList(txID)
        % Normalize channel
        normChan = normalize(chTot(txID,:));
        % Apply interested channel samples
        channel = normChan(:,tol:maxIter-1);  % Take one sample at the end
        channel_aps = ones(size(channel,maxIter-tol));
        channel_phs = angle(channel);
        % Retrieve the applied channel
        channel_est_aps = ones(size(channel,maxIter-tol));
        channel_est_phs = angle(normChan(:,tol+1:maxIter)) + pi;
        % Aggregate the capacities
        ch1 = ch1 + abs(cos(channel_phs - channel_phs));
        ch2 = ch2 + abs(cos(channel_phs + channel_est_phs));
    end    
    cap_max(:,txID) = log(1 + power_lin.*ch1);
    cap_achieved(:,txID) = log(1 + power_lin.*ch2);
%     subplot(length(nTxAntennasList),1,txID); hold on;
%     plot((1:size(cap_max,1)),cap_max(:,txID),'LineWidth',2,'Color','b');
%     plot((1:size(cap_max,1)),cap_achieved(:,txID),'LineWidth',2,'Color','r');
%     ylim([0 3]);
%     xlabel('Iterations','FontSize',12);
%     ylabel('Mutual info','FontSize',12);
%     title(strcat('NTx:',{' '},num2str(nTxAnt)),'FontSize',12);
    
    axes(ha(txID)); hold on;
    plot((1:size(cap_max,1)),cap_max(:,txID),'LineWidth',2,'Marker','none','Color','b');
    plot((1:size(cap_max,1)),cap_achieved(:,txID),'LineWidth',2,'Marker','none','Color','r');
    xlab = xlabel('Iteration','FontSize',12);
    % Manipulate y-axis
    yyaxis left
    ylabel(strcat('NTx:',{' '},num2str(nTxAnt)),'FontSize',12,'FontWeight','bold','Color','k');
    ylim([0 3]);
    yyaxis right
    ylab = ylabel({'Mutual';'Information'},'FontSize',12,'FontWeight','bold','Color','k');
    set(ylab,'Rotation',0,'VerticalAlignment','middle', 'HorizontalAlignment','left');
%     set(gca,'ytick',[]);
%     p = get(xlab,'position'); p(2) = 0.9*p(2); set(xlab,'position',p);
%     p1 = p.*[2.15 .5 1];  % Horizontal x Vertical x Depth
%     set(xlab,'Position',p1);
%     set(xlab);
%     set(ha(txID),'XTickLabel',''); set(ha,'YTickLabel','')
    grid minor;
end
leg = legend('Ideal','Applied');
set(leg,'FontSize',12,'Orientation','Horizontal');

%%
load('BERTot_file','BER');
% load('BERTot_Rx_file','BER_RxHov');
% BER = BER_RxHov;
% load('BERTot_Tx_KRI_file','BER_TxHov_KRI');
% BER = BER_TxHov_KRI;

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
title('BER improvement with Beamforming - UAVs flying','FontSize',14);