function BER_receiverPlot(varargin)

addpath('kriging/');  % Include Kriging folder (variogram and kriging)
figIdx = 20;

if nargin==0
%     results = 'sim_BER-exp2ant.mat';
%     results = 'res_ROT_indoors.mat';
    results = 'results/res_ROT_outdoors.mat';
elseif nargin==1
    results = varargin{1};
else
    error('Wrong input parameters');
end
    
load(results);  % Load results
% lastIter = i;

%% Plot Channel Gain in 2D
% figure(figIdx);  figIdx = figIdx + 1;
% for id = 1:nTxAntennas
%     subplot(1,nTxAntennas,id); hold on;
%     yyaxis left
%     plot((1:lastIter),abs(chTot(id,1:lastIter)),'Color','k','LineWidth',2.5,'Marker','none','MarkerSize',2);
%     ylabel('Gain','FontSize',12);
%     yyaxis right
%     plot((1:lastIter),appliedElev(1:lastIter),'Color','r','LineWidth',1.5,'LineStyle','-');
%     plot((1:lastIter),appliedAzym(1:lastIter),'Color','b','LineWidth',1.5,'LineStyle','-');
%     title('Absolute Gain','FontSize',12);
%     xlabel('Iteration','FontSize',12);
%     ylabel('Angle (degree)','FontSize',12);
%     lg = legend('Gain','Elevation angle','Azimuth angle');
%     set(lg,'FontSize',12);
% end

%% Plot Channel Gain in 3D
chGain = abs(chTot);
[X,Y] = meshgrid(elevList,azimList);
dim = size(X,1)*size(X,2);

final_av = zeros(nTxAntennas,length(elevList)*length(azimList));
final_var = zeros(nTxAntennas,length(elevList)*length(azimList));
final_av_tot = zeros(length(azimList),length(elevList),nTxAntennas);
final_var_tot = zeros(length(azimList),length(elevList),nTxAntennas);
leg = cell(nTxAntennas,1);
for id = 1:nTxAntennas
    for t = 1:dim
        elev = X(t);
        azym = Y(t);
        idx_elev = find(appliedElev==elev);
        idx_azym = find(appliedAzym==azym);
        idx = intersect(idx_elev,idx_azym);
        final_av(id,t) = mean(chGain(id,idx));
        final_var(id,t) = var(chGain(id,idx));
    end
    final_av_tot(:,:,id) = reshape(final_av(id,:),[size(X,1),size(X,2)]);
    final_var_tot(:,:,id) = reshape(final_var(id,:),[size(X,1),size(X,2)]);
    elevList_int = (max(elevList):-5:min(elevList));
    azimList_int = (max(azimList):-5:min(azimList));
    [X1,Y1] = meshgrid(elevList_int,azimList_int);
    Z1 = interp2(X,Y,final_av_tot(:,:,id),X1,Y1,'spline');
    % Plot Empirical Gain
    figure(figIdx + id);
    surf(X1,Y1,Z1);
    xlabel('Elevation','FontSize',12);
    ylabel('Azimuth','FontSize',12);
    title(strcat('antenna #',{' '},num2str(id)),'FontSize',12);
    % Plot Statistics
    A = reshape(X1,size(X1,1)*size(X1,2),1);
    B = reshape(Y1,size(Y1,1)*size(Y1,2),1);
    C = reshape(Z1,size(Z1,1)*size(Z1,2),1);
    d2 = variogram([A B],normalize(C),'plotit',false,'nrbins',30);
    figure(figIdx + nTxAntennas + 1);  hold on;
    plot(d2.distance,d2.val,'LineWidth',1.5,'LineStyle','--','Marker','o');
    xlabel('h','FontSize',12);
    ylabel('\gamma (h)','FontSize',12);
    title(strcat('(Semi-)Variogram antenna #',{' '},num2str(id)),'FontSize',12);
    leg(id) = strcat('antenna #',{' '},num2str(id));
    figure(figIdx + nTxAntennas + 2);  hold on;
    [~,~,~,~] = variogramfit(d2.distance,d2.val,[],[],[],'model','gaussian');
end
figure(figIdx + nTxAntennas + 1);
legend(leg,'FontSize',12,'Location','NorthWest');
figIdx = figIdx + 1;

%% Plot Gain (Real, Imag and Absolute) for each antenna
% for id = 1:nTxAntennas
%     figure(figIdx);
%     subplot(nTxAntennas,3,3*(id-1) + 1); hold on; grid minor;
%     % plot((1:maxIter),real(chRealTot),'Color','b');
%     plot((1:lastIter),real(chTot(id,1:lastIter)),'Color','r','LineWidth',2);
%     title('Real','FontSize',12);
%     xlabel('Iteration','FontSize',12);
%     ylabel('Gain','FontSize',12);
%     subplot(nTxAntennas,3,3*(id-1) + 2); hold on; grid minor;
%     % plot((1:maxIter),imag(chRealTot),'Color','b');
%     plot((1:lastIter),imag(chTot(id,1:lastIter)),'Color','r','LineWidth',2);
%     title('Imaginary','FontSize',12);
%     xlabel('Iteration','FontSize',12);
%     ylabel('Gain','FontSize',12);
%     subplot(nTxAntennas,3,3*(id-1) + 3); hold on; grid minor;
%     plot((1:lastIter),abs(chTot(id,1:lastIter)),'Color','r','LineWidth',2);
%     title('Absolute Gain','FontSize',12);
%     xlabel('Iteration','FontSize',12);
%     ylabel('Gain','FontSize',12);
% 
%     figure(figIdx + 1);
%     subplot(4,1,id); hold on;
%     covTot_imag = zeros(lastIter,1);
%     covTot_real = zeros(lastIter,1);
%     re = real(chTot(id,:));
%     im = imag(chTot(id,:));
%     re(isnan(re))=0;
%     window_size = 1;
%     for k = 1:1:lastIter-window_size
%         covTot_imag(k) = cov(im(k:k+window_size));
%         covTot_real(k) = cov(re(k:k+window_size));
%     end
%     plot(abs(1-covTot_imag),'LineWidth',2);
%     plot(abs(1-covTot_real),'LineWidth',2);
%     legend('Imaginary','Real');
%     ylim([0.999 1]);
% end
% figIdx = figIdx + 1;

%% Plot cumulative Gain across antennas
% figure(figIdx);  figIdx = figIdx + 1;
% totGain = 0;
% for txID = 1:nTxAntennas
%     t = abs(chTot(txID,1:lastIter));
%     t(isnan(t))=0;  % replace nan values
%     totGain = totGain + t;
% end
%  plot((1:lastIter),totGain,'Color','r','LineWidth',2);
% title('Total channel Gain','FontSize',12);
% xlabel('Iteration','FontSize',12);
% ylabel('Gain (linear)','FontSize',12);

%% Plot measured BER
% figure(figIdx);  figIdx = figIdx + 1;
% subplot(611); plot(abs(BER(1:lastIter,1)),'LineWidth',1.5); legend('64-QAM'); ylabel('BER'); ylim([0 1]);
% subplot(612); plot(abs(BER(1:lastIter,2)),'LineWidth',1.5); legend('32-QAM'); ylabel('BER'); ylim([0 1]);
% subplot(613); plot(abs(BER(1:lastIter,3)),'LineWidth',1.5); legend('16-QAM'); ylabel('BER'); ylim([0 1]);
% subplot(614); plot(abs(BER(1:lastIter,4)),'LineWidth',1.5); legend('8-QAM'); ylabel('BER'); ylim([0 1]);
% subplot(615); plot(abs(BER(1:lastIter,5)),'LineWidth',1.5); legend('QPSK'); ylabel('BER'); ylim([0 1]);
% subplot(616); plot(abs(BER(1:lastIter,6)),'LineWidth',1.5); legend('BPSK'); ylabel('BER'); ylim([0 1]);
