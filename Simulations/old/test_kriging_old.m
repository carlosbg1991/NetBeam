clear all; clear classes;  %#ok
% close all; clc;
addpath('kriging/');  % Include Kriging folder (variogram and kriging)
addpath('BrewerMap/');  % Include additional colors: 'BrBG'|'PRGn'|'PiYG'|'PuOr'|'RdBu'|'RdGy'|'RdYlBu'|'RdYlGn'|'Spectral'

% load('results/res_ROT_outdoors.mat','azimList','elevList','chTot','nTxAntennas','appliedElev','appliedAzym');
% load('results-19-Nov-2018_16:17:56.mat','azimList','elevList','chTot','nTxAntennas','appliedElev','appliedAzym');
% load('results-19-Nov-2018_18:52:28.mat','azimList','elevList','chTot','nTxAntennas','appliedElev','appliedAzym');
% load('results-19-Nov-2018_18:52:28_2.mat','azimList','elevList','chTot','nTxAntennas','appliedElev','appliedAzym');
% load('results-25-Nov-2018_14:49:13.mat','azimList','elevList','chTot','nTxAntennas','appliedElev','appliedAzym');
% load('results-25-Nov-2018_14:55:08.mat','azimList','elevList','chTot','nTxAntennas','appliedElev','appliedAzym');
% load('results-25-Nov-2018_20:15:13.mat','azimList','elevList','chTot','nTxAntennas','appliedElev','appliedAzym');

load('results-27-Nov-2018_01:32:43.mat','azimList','elevList','chTot','nTxAntennas','appliedElev','appliedAzym');

% load('results-29-Nov-2018_04:41:46.mat','azimList','elevList','chTot','nTxAntennas','appliedElev','appliedAzym');
% load('results-29-Nov-2018_22:45:13.mat','azimList','elevList','chTot','nTxAntennas','appliedElev','appliedAzym');
% load('results-03-Dec-2018_13:12:23.mat','azimList','elevList','chTot','nTxAntennas','appliedElev','appliedAzym');
% load('results-03-Dec-2018_14:22:34.mat','azimList','elevList','chTot','nTxAntennas','appliedElev','appliedAzym');
% load('results-03-Dec-2018_15:32:52.mat','azimList','elevList','chTot','nTxAntennas','appliedElev','appliedAzym');
% load('results-03-Dec-2018_16:41:21.mat','azimList','elevList','chTot','nTxAntennas','appliedElev','appliedAzym');
% load('results-03-Dec-2018_17:49:02.mat','azimList','elevList','chTot','nTxAntennas','appliedElev','appliedAzym');
% load('results-04-Dec-2018_15:05:21.mat','azimList','elevList','chTot','nTxAntennas','appliedElev','appliedAzym');

%% Parse inputs
chGain = abs(chTot);
[X,Y] = meshgrid(elevList,azimList);
dim = size(X,1)*size(X,2);
plotIdx = 1;

%% Parse exhaustive
Z_prel = zeros(nTxAntennas,length(elevList)*length(azimList));
Z = zeros(length(azimList),length(elevList),nTxAntennas);
for id = 1:nTxAntennas
    for t = 1:dim
        elev = X(t);
        azym = Y(t);
        idx_elev = find(appliedElev==elev);
        idx_azym = find(appliedAzym==azym);
        idx = intersect(idx_elev,idx_azym);
        Z_prel(id,t) = mean(chGain(id,idx));
    end
    Z(:,:,id) = reshape(Z_prel(id,:),[size(X,1),size(X,2)]);
end

% %% Test Fixed 
% n_samples = 10;
% for id = 1:nTxAntennas
%     % New interpolated input (assume this is the real exhaustive)
%     elevList_int = (max(elevList):-1:min(elevList));
%     azimList_int = (max(azimList):-1:min(azimList));
%     [X0,Y0] = meshgrid(elevList_int,azimList_int);
%     Z0(:,:,id) = interp2(X,Y,Z(:,:,id),X0,Y0,'spline');  %#ok
% 
%     % Generate Exhaustive image
%     subplot(nTxAntennas,4,plotIdx); hold on;  plotIdx = plotIdx + 1;
%     imagesc(X0(1,:),Y0(:,1),Z0(:,:,id)); axis image; axis xy
% 
%     % Sample the space
%     combo = combvec(elevList_int,azimList_int);  % Truth table
%     samp_ids = randperm(length(combo));  % Generate Random sequence (no EI for now)
%     samples = combo(:,samp_ids(1:n_samples));  % Extract samples Elevations and Azimuths
%     elev_sample = samples(1,:);
%     azim_sample = samples(2,:);
%     
% %     % Combo 1    
% %     elev_sample = [20 20 20 60 60 60 80 80 80];
% %     azim_sample = [40 90 140 40 90 140 40 90 140];
% %     % Combo 2
% %     elev_sample = [20 20 80 80];
% %     azim_sample = [40 140 40 140];
%     
%     [X_sample,Y_sample] = meshgrid(elev_sample,azim_sample);
%     Z_sample = interp2(X0,Y0,Z0(:,:,id),X_sample,Y_sample);
%     plot(elev_sample,azim_sample,'+r','MarkerSize',5);
%     title('Real Map and sampling');
% 
%     % calculate the sample variogram
%     a = reshape(X_sample,[size(X_sample,1)*size(X_sample,2),1]);
%     b = reshape(Y_sample,[size(Y_sample,1)*size(Y_sample,2),1]);
%     c = reshape(Z_sample,[size(Z_sample,1)*size(Z_sample,2),1]);
%     v = variogram([a b],c,'plotit',false,'maxdist',100);
% 
%     % Compute Gaussian fit following Least Squares
%      subplot(nTxAntennas,4,plotIdx); hold on;  plotIdx = plotIdx + 1;
%     [~,~,~,vstruct] = variogramfit(v.distance,v.val,[],[],[],'model','gaussian');
%     title('Prior fitting: Spatial correlation');
% 
%     % now use the sampled locations in a kriging
%     subplot(nTxAntennas,4,plotIdx); hold on;  plotIdx = plotIdx + 1;
%     [Zhat,Zvar] = kriging(vstruct,X_sample,Y_sample,Z_sample,X0,Y0);
%     imagesc(X0(1,:),Y0(:,1),Zhat); axis image; axis xy
%     title('kriging predictions')
% 
%     subplot(nTxAntennas,4,plotIdx); hold on;  plotIdx = plotIdx + 1;
%     contour(X0,Y0,Zvar); axis image
% %     contourf(Zvar);
%     title('kriging variance')
% end

%% Test Progressive
n_samples     = 50;    % Maximum number of samples to run Kriging with
minSamples    = 15;     % Minimum Kriging initialization space
new_additions = 1;     % New trials per iteration
newAddPoss    = 4;     % New possible possitions for DIRECT
id            = 3;     % For instance, could be 1,2,3,4 (antenna id)
plotFlag      = true;  % Flag to plot results

nSamplesList = minSamples:new_additions:n_samples;  % List of samples to try

% New interpolated input (assume this is the real exhaustive)
elevList_int = (max(elevList):-1:min(elevList));
azimList_int = (max(azimList):-1:min(azimList));
[X0,Y0] = meshgrid(elevList_int,azimList_int);
Z0 = interp2(X,Y,Z(:,:,id),X0,Y0,'spline');

% Performance indicators - VARIANCE
avVar_random = zeros(size(nSamplesList));
maxVar_random = zeros(size(nSamplesList));
minVar_random = zeros(size(nSamplesList));
avVar_DIRECTrandom = zeros(size(nSamplesList));
maxVar_DIRECTrandom = zeros(size(nSamplesList));
minVar_DIRECTrandom = zeros(size(nSamplesList));
avVar_DIRECTminVar = zeros(size(nSamplesList));
maxVar_DIRECTminVar = zeros(size(nSamplesList));
minVar_DIRECTminVar = zeros(size(nSamplesList));

% Performance indicators - OPTIMUM ANGLE PREDICTION and GAIN PREDICTION
pred_azim_random = zeros(size(nSamplesList));
pred_azim_DIRECTrandom = zeros(size(nSamplesList));
pred_azim_DIRECTminVar = zeros(size(nSamplesList));
pred_elev_random = zeros(size(nSamplesList));
pred_elev_DIRECTrandom = zeros(size(nSamplesList));
pred_elev_DIRECTminVar = zeros(size(nSamplesList));
pred_gain_random = zeros(size(nSamplesList));
pred_gain_DIRECTrandom = zeros(size(nSamplesList));
pred_gain_DIRECTminVar = zeros(size(nSamplesList));

% Performance indicators - RELIABILITY on predicted MAP
missmatch_random = zeros(size(nSamplesList));
missmatch_DIRECTrandom = zeros(size(nSamplesList));
missmatch_DIRECTminVar = zeros(size(nSamplesList));

% find optimum angles - Exhaustive
[exhv_gain, max_idx] = max(Z0(:));
[indexX,indexY] = ind2sub(size(Z0),max_idx);
exhv_elevation = X0(indexX,indexY);
exhv_azimuth = Y0(indexX,indexY);
pred_azim_exhaust = repmat(exhv_azimuth,size(nSamplesList));
pred_elev_exhaust = repmat(exhv_elevation,size(nSamplesList));
pred_gain_exhaust = repmat(exhv_gain,size(nSamplesList));

selPolicyList = {'DIRECT-minVar','random','DIRECT-rand'}; % Chose between 'random', 'EI', 'PI', 'DIRECT-rand', 'DIRECT-minVar'
for pol = selPolicyList
    
    fprintf('Executing with policy "%s"\n',pol{:});
    selPolicy = pol{:};
    
    % Initialize Sampling
    if strcmp(selPolicy,'random')
        % Initialize random space sampling
        combo = combvec(elevList_int,azimList_int);  % Truth table
        samp_ids = randperm(length(combo));  % Generate Random sequence (no EI for now)
        samples = combo(:,samp_ids(1:n_samples));  % Extract samples Elevations and Azimuths
        elev_sample_random = samples(1,:);
        azim_sample_random = samples(2,:);
        elev_sample1 = elev_sample_random(1:minSamples);  % Store the sampled elevation angles
        azim_sample1 = azim_sample_random(1:minSamples);  % Store the sampled azimuth angles
    elseif strcmp(selPolicy,'EI')
        % To-do
    elseif strcmp(selPolicy,'PI')
        % To-do
    elseif strcmp(selPolicy,'DIRECT-rand') || strcmp(selPolicy,'DIRECT-minVar')
        % Initialize DIRECT
        selPoss = CBG_DIRECT([0 180], [0 90], minSamples);
        % Plot Initial DIRECT configuration
        if plotFlag
            figure(2); cla reset; hold on;
            scatter(selPoss(1,1,:),selPoss(1,2,:))
            for t = 1:newAddPoss
                plot([selPoss(2,1,t) selPoss(2,1,t)],[selPoss(3,1,t) selPoss(3,2,t)],'lineWidth',2)
                plot([selPoss(2,1,t) selPoss(2,2,t)],[selPoss(3,1,t) selPoss(3,1,t)],'lineWidth',2)
            end
            hold off;
        end
        % First, select the samples to trye (special case, more than one)
        elev_sample1 = selPoss(1,2,:);  elev_sample1 = elev_sample1(:).';
        azim_sample1 = selPoss(1,1,:);  azim_sample1 = azim_sample1(:).';
        % Second, initialize massively DIRECT (special case, more than one)
        for selTrial = 1:newAddPoss
            selCenter = selPoss(1,:,selTrial);
            selLimAzim = selPoss(2,:,selTrial);
            selLimElev = selPoss(3,:,selTrial);
            % Second, break down the new area into new sub-areas using DIRECT
            addSelPoss = CBG_DIRECT(selLimAzim,selLimElev,newAddPoss);
            % Append new sub-area to global area
            selPoss = cat(3, selPoss, addSelPoss);
        end
        % Remove explored possibility from selection set
        selPoss(:,:,1:newAddPoss) = [];
    end

    count = 1;  % To index results into storing
    for iter = nSamplesList

        [X_sample,Y_sample] = meshgrid(elev_sample1,azim_sample1);
        Z_sample = interp2(X0,Y0,Z0,X_sample,Y_sample);

        % calculate the sample variogram
        a = reshape(X_sample,[size(X_sample,1)*size(X_sample,2),1]);
        b = reshape(Y_sample,[size(Y_sample,1)*size(Y_sample,2),1]);
        c = reshape(Z_sample,[size(Z_sample,1)*size(Z_sample,2),1]);
        v = variogram([a b],c,'maxdist',100,'plotit',false);

        % compute Gaussian fit following Least Squares
        [~,~,~,vstruct] = variogramfit(v.distance,v.val,[],[],[],'model','gaussian','plotit',false);

        % use the sampled locations in a kriging
        [Zhat,Zvar] = kriging(vstruct,X_sample,Y_sample,Z_sample,X0,Y0);
        
        % find optimum angles - prediction
        [pred_gain, max_idx] = max(Zhat(:));
        [X,Y] = ind2sub(size(Zhat),max_idx);
        pred_elevation = X0(X,Y);
        pred_azimuth = Y0(X,Y);

        if plotFlag
            
            % Plot DIRECT map
            if strcmp(selPolicy,'DIRECT-rand') || strcmp(selPolicy,'DIRECT-minVar')
                % Plot possibilities
                figure(2); cla reset; hold on;
                scatter(selPoss(1,1,:),selPoss(1,2,:))
                for t = 1:size(selPoss,3)
                    plot([selPoss(2,1,t) selPoss(2,1,t)],[selPoss(3,1,t) selPoss(3,2,t)],'lineWidth',2)
                    plot([selPoss(2,1,t) selPoss(2,2,t)],[selPoss(3,1,t) selPoss(3,1,t)],'lineWidth',2)
                end
            end

            % Plot exhaustive, predictions and uncertainty
            figure(1); 

            ax = subplot(1,4,1); cla reset; hold on;
            brewermap('*RdBu');
            colormap(ax,'brewermap');
            colorbar('eastoutside');
            imagesc(X0(1,:),Y0(:,1),Z0); axis image; axis xy  % Generate Exhaustive image
            plot(elev_sample1,azim_sample1,'+r','MarkerSize',5);  % Plot trials so far
            mystr = strcat('#',num2str(iter),{' '},'Real+sampling');
            title(mystr);
            hold off

            ax = subplot(1,4,2); cla reset;
            imagesc(X0(1,:),Y0(:,1),Zhat); axis image; axis xy;
            brewermap('*RdBu');
            colormap(ax,'brewermap');
            colorbar('eastoutside');
            mystr = strcat('#',num2str(iter),{' '},'kriging predictions');
            title(mystr);

            ax = subplot(1,4,3); cla reset;
            contourf(X0,Y0,Zvar);
            mystr = strcat('#',num2str(iter),{' '},'kriging variance');
            title(mystr);
            brewermap('*RdBu');
            colormap(ax,'brewermap');
            colorbar('eastoutside');

            ax = subplot(1,4,4); cla reset;
            plot(abs(Zvar(:)));
            brewermap('*RdBu');
            colormap(ax,'brewermap');
            colorbar('eastoutside');
            mystr = strcat('#',num2str(iter),{' '},'Overall variance');
            title(mystr);
        end

        % Select next trial upong selection policy
        if strcmp(selPolicy,'random')
            elev_sample1 = elev_sample_random(1:iter);
            azim_sample1 = azim_sample_random(1:iter);
            % Store results - variance
            avVar_random(count) = mean(abs(Zvar(:)));
            maxVar_random(count) = max(abs(Zvar(:)));
            minVar_random(count) = min(abs(Zvar(:)));
            % Store results - predictions
            pred_azim_random(count) = pred_azimuth;
            pred_elev_random(count) = pred_elevation;
            pred_gain_random(count) = pred_gain;
        elseif strcmp(selPolicy,'EI')
            % To-do
        elseif strcmp(selPolicy,'PI')
            % To-do
        elseif strcmp(selPolicy,'DIRECT-rand')
            % First, select sample following a random selection
            possIndices = size(selPoss,3);
            selTrial = randi([1 possIndices]);
            % Retrieve results
            selCenter = selPoss(1,:,selTrial);
            selLimAzim = selPoss(2,:,selTrial);
            selLimElev = selPoss(3,:,selTrial);
            % Second, break down the new area into new sub-areas using DIRECT
            addSelPoss = CBG_DIRECT(selLimAzim,selLimElev,newAddPoss);
            % Remove explored possibility from selection set
            selPoss(:,:,selTrial) = [];
            % Append new sub-area to global area
            selPoss = cat(3, selPoss, addSelPoss);
            % Select actual azimuth and elevation to evaluate with kriging
            elev_sample1 = [elev_sample1 selCenter(2)];  %#ok<AGROW>
            azim_sample1 = [azim_sample1 selCenter(1)];  %#ok<AGROW>
            % Store results - variance
            avVar_DIRECTrandom(count) = mean(abs(Zvar(:)));
            maxVar_DIRECTrandom(count) = max(abs(Zvar(:)));
            minVar_DIRECTrandom(count) = min(abs(Zvar(:)));
            % Store results - predictions
            pred_azim_DIRECTrandom(count) = pred_azimuth;
            pred_elev_DIRECTrandom(count) = pred_elevation;
            pred_gain_DIRECTrandom(count) = pred_gain;
        elseif strcmp(selPolicy,'DIRECT-minVar')
            % First, select sample with lowest variance
            myAzims = selPoss(1,1,:);   myAzims = myAzims(:);
            myElevs = selPoss(1,2,:);   myElevs = myElevs(:);
            myX0 = X0(:);  myY0 = Y0(:);  myZvar = Zvar(:);
            myY01 = repmat(myY0,1,length(myAzims)).';
            myX01 = repmat(myX0,1,length(myElevs)).';
            myIndices = zeros(1,length(myAzims));
            for temp = 1:length(myAzims)
                test_azim = find(myY0==ceil(myAzims(temp)));  %Swap them, needed
                test_elev = find(myX0==ceil(myElevs(temp)));  %Swap them, needed
                myIndices(temp) = intersect(test_azim,test_elev);
            end
            [~,selTrial] = max(myZvar(myIndices));
            % Retrieve results
            selCenter = selPoss(1,:,selTrial);
            selLimAzim = selPoss(2,:,selTrial);
            selLimElev = selPoss(3,:,selTrial);
            % Second, break down the new area into new sub-areas using DIRECT
            addSelPoss = CBG_DIRECT(selLimAzim,selLimElev,newAddPoss);
            % Remove explored possibility from selection set
            selPoss(:,:,selTrial) = [];
            % Append new sub-area to global area
            selPoss = cat(3, selPoss, addSelPoss);
            % Select actual azimuth and elevation to evaluate with kriging
            elev_sample1 = [elev_sample1 selCenter(2)];  %#ok<AGROW>
            azim_sample1 = [azim_sample1 selCenter(1)];  %#ok<AGROW>
            % Store results - variance
            avVar_DIRECTminVar(count) = mean(abs(Zvar(:)));
            maxVar_DIRECTminVar(count) = max(abs(Zvar(:)));
            minVar_DIRECTminVar(count) = min(abs(Zvar(:)));
            % Store results - predictions
            pred_azim_DIRECTminVar(count) = pred_azimuth;
            pred_elev_DIRECTminVar(count) = pred_elevation;
            pred_gain_DIRECTminVar(count) = pred_gain;
        end
        
        % Update counter for results
        count = count + 1;
    end
    % Clean intermediate variables
end

figure(3); hold on;
plot(nSamplesList,avVar_random,'lineWidth',2);
plot(nSamplesList,avVar_DIRECTrandom,'lineWidth',2);
plot(nSamplesList,avVar_DIRECTminVar,'lineWidth',2);
grid minor;
lg = legend('RANDOM-samp and RANDOM-sel','DIRECT-samp and RANDOM-sel','DIRECT-samp and MINVAR-sel');
set(lg,'FontSize',12);
title('analysis on the AVERAGE Variance achieved','FontSize',12);
xlabel('Number of trials','FontSize',12);
ylabel('Variance','FontSize',12);
hold off;

figure(4); hold on;
plot(nSamplesList,maxVar_random,'lineWidth',2);
plot(nSamplesList,maxVar_DIRECTrandom,'lineWidth',2);
plot(nSamplesList,maxVar_DIRECTminVar,'lineWidth',2);
grid minor;
lg = legend('RANDOM-samp and RANDOM-sel','DIRECT-samp and RANDOM-sel','DIRECT-samp and MINVAR-sel');
set(lg,'FontSize',12);
title('analysis on the MAX Variance achieved','FontSize',12);
xlabel('Number of trials','FontSize',12);
ylabel('Variance','FontSize',12);

figure(5); hold on;
plot(nSamplesList,minVar_random,'lineWidth',2);
plot(nSamplesList,minVar_DIRECTrandom,'lineWidth',2);
plot(nSamplesList,minVar_DIRECTminVar,'lineWidth',2);
grid minor;
lg = legend('RANDOM-samp and RANDOM-sel','DIRECT-samp and RANDOM-sel','DIRECT-samp and MINVAR-sel');
set(lg,'FontSize',12);
title('analysis on the MIN Variance achieved','FontSize',12);
xlabel('Number of trials','FontSize',12);
ylabel('Variance','FontSize',12);
hold off

figure(6); hold on;
plot(nSamplesList,pred_gain_random,'lineWidth',2);
plot(nSamplesList,pred_gain_DIRECTrandom,'lineWidth',2);
plot(nSamplesList,pred_gain_DIRECTminVar,'lineWidth',2);
plot(nSamplesList,exhv_gain,'lineWidth',2);
grid minor;
lg = legend('RANDOM-samp and RANDOM-sel','DIRECT-samp and RANDOM-sel','DIRECT-samp and MINVAR-sel','Optimum Exhaustive');
set(lg,'FontSize',12);
title('analysis on the MAX GAIN achieved','FontSize',12);
xlabel('Number of trials','FontSize',12);
ylabel('Variance','FontSize',12);
hold off;

figure(7); hold on;
plot(nSamplesList,pred_azim_random,'lineWidth',2);
plot(nSamplesList,pred_azim_DIRECTrandom,'lineWidth',2);
plot(nSamplesList,pred_azim_DIRECTminVar,'lineWidth',2);
plot(nSamplesList,exhv_azimuth,'lineWidth',2);
grid minor;
lg = legend('RANDOM-samp and RANDOM-sel','DIRECT-samp and RANDOM-sel','DIRECT-samp and MINVAR-sel','Optimum Exhaustive');
set(lg,'FontSize',12);
title('analysis on the AZIMUTH prediction','FontSize',12);
xlabel('Number of trials','FontSize',12);
ylabel('Variance','FontSize',12);
hold off;

figure(8); hold on;
plot(nSamplesList,pred_elev_random,'lineWidth',2);
plot(nSamplesList,pred_elev_DIRECTrandom,'lineWidth',2);
plot(nSamplesList,pred_elev_DIRECTminVar,'lineWidth',2);
plot(nSamplesList,exhv_elevation,'lineWidth',2);
grid minor;
lg = legend('RANDOM-samp and RANDOM-sel','DIRECT-samp and RANDOM-sel','DIRECT-samp and MINVAR-sel','Optimum Exhaustive');
set(lg,'FontSize',12);
title('analysis on the AZIMUTH prediction','FontSize',12);
xlabel('Number of trials','FontSize',12);
ylabel('Variance','FontSize',12);
hold off;