% clear all; clear classes;  %#ok
close all; clc;
addpath('kriging/');  % Include Kriging folder (variogram and kriging)
addpath('BrewerMap/');  % Include additional colors: 'BrBG'|'PRGn'|'PiYG'|'PuOr'|'RdBu'|'RdGy'|'RdYlBu'|'RdYlGn'|'Spectral'
set(0,'DefaultFigureColor','remove');  % No gray background in figures

%% PARAMETERS
minSamples    = 4;     % Minimum Kriging initialization space
n_samples     = 8;    % Maximum number of samples to run Kriging with
new_additions = 1;     % New trials per iteration
newAddPoss    = 4;     % New possible possitions for DIRECT
nIter         = 1;    % Iterations to run Random over
antID         = 2;     % Antenna ID, could be 1,2,3,4
expID         = 2;     % Experiment ID, could be 1,2,3,4,5
plotFlag      = true;  % Flag to plot results

%% PARSE Data if not done before
if ~exist('outdoor','var') || ~exist('indoor','var')
    CBG_parse_experiments;  % Parse data into structs 'indoor', 'outdoor'
end
elevList = indoor.elevList;
azimList = indoor.azimList;
[X,Y] = meshgrid(elevList,azimList);  % Orientation space (exhaustive)
Z = indoor.gainTot(:,:,antID,expID);  % Channel gain (exhaustive)

nSamplesList = minSamples:new_additions:n_samples;  % List of samples to try

% New interpolated input (assume this is the real exhaustive)
elevList_int = (max(elevList):-1:min(elevList));
azimList_int = (max(azimList):-1:min(azimList));
[X0,Y0] = meshgrid(elevList_int,azimList_int);
Z0 = interp2(X,Y,Z,X0,Y0,'spline');

% Performance indicators - VARIANCE
avVar_random = zeros(size(nSamplesList));  % DIRECT random
maxVar_random = zeros(size(nSamplesList));
minVar_random = zeros(size(nSamplesList));
avVar_DIRECTrandom = zeros(size(nSamplesList));  % DIRECT random  
maxVar_DIRECTrandom = zeros(size(nSamplesList));
minVar_DIRECTrandom = zeros(size(nSamplesList));
avVar_DIRECTminVar = zeros(size(nSamplesList));  % DIRECT minVar  
maxVar_DIRECTminVar = zeros(size(nSamplesList));
minVar_DIRECTminVar = zeros(size(nSamplesList));
avVar_UM = zeros(size(nSamplesList));  % Uncertainty minimization
maxVar_UM = zeros(size(nSamplesList));
minVar_UM = zeros(size(nSamplesList));

% Performance indicators - OPTIMUM ANGLE PREDICTION and GAIN PREDICTION
ml_random_pred_azim = zeros(size(nSamplesList));  % Random
ml_random_pred_elev = zeros(size(nSamplesList));
ml_random_pred_gain = zeros(size(nSamplesList));
ml_random_real_gain = zeros(size(nSamplesList));
ml_random_pred_gap_gain = zeros(size(nSamplesList));
ml_random_real_gap_gain = zeros(size(nSamplesList));
ml_DIRECTrandom_pred_azim = zeros(size(nSamplesList));  % DIRECT random  
ml_DIRECTrandom_pred_elev = zeros(size(nSamplesList));
ml_DIRECTrandom_pred_gain = zeros(size(nSamplesList));
ml_DIRECTrandom_real_gain = zeros(size(nSamplesList));
ml_DIRECTrandom_pred_gap_gain = zeros(size(nSamplesList));
ml_DIRECTrandom_real_gap_gain = zeros(size(nSamplesList));
ml_DIRECTminVar_pred_azim = zeros(size(nSamplesList));  % DIRECT minVar  
ml_DIRECTminVar_pred_elev = zeros(size(nSamplesList));
ml_DIRECTminVar_pred_gain = zeros(size(nSamplesList));
ml_DIRECTminVar_real_gain = zeros(size(nSamplesList));
ml_DIRECTminVar_pred_gap_gain = zeros(size(nSamplesList));
ml_DIRECTminVar_real_gap_gain = zeros(size(nSamplesList));
ml_UM_pred_azim = zeros(size(nSamplesList));  % Uncertainty minimization
ml_UM_pred_elev = zeros(size(nSamplesList));
ml_UM_pred_gain = zeros(size(nSamplesList));
ml_UM_real_gain = zeros(size(nSamplesList));
ml_UM_pred_gap_gain = zeros(size(nSamplesList));
ml_UM_real_gap_gain = zeros(size(nSamplesList));

% Extra results for Random 
avVar_random1 = zeros(length(nSamplesList),nIter);
maxVar_random1 = zeros(length(nSamplesList),nIter);
minVar_random1 = zeros(length(nSamplesList),nIter);
ml_random_pred_gain1 = zeros(length(nSamplesList),nIter);
ml_random_real_gain1 = zeros(length(nSamplesList),nIter);

% Performance indicators - RELIABILITY on predicted MAP
missmatch_random = zeros(size(nSamplesList));
missmatch_DIRECTrandom = zeros(size(nSamplesList));
missmatch_DIRECTminVar = zeros(size(nSamplesList));

% find optimum angles - Exhaustive
[exhv_gain, max_idx] = max(Z0(:));
[idx_exhv_maxGain_X,idx_exhv_maxGain_Y] = ind2sub(size(Z0),max_idx);
exhv_elevation = X0(idx_exhv_maxGain_X,idx_exhv_maxGain_Y);
exhv_azimuth = Y0(idx_exhv_maxGain_X,idx_exhv_maxGain_Y);
exhv_azimtuh_long = repmat(exhv_azimuth,size(nSamplesList));
exhv_elevation_long = repmat(exhv_elevation,size(nSamplesList));
exhv_gain_long = repmat(exhv_gain,size(nSamplesList));

% selPolicyList = {'DIRECT-minVar','random','UM','DIRECT-rand'}; % Chose between 'random', 'EI', 'PI', 'UM', 'DIRECT-rand', 'DIRECT-minVar'
selPolicyList = {'PI'}; % Chose between 'random', 'EI', 'PI', 'UM', 'DIRECT-rand', 'DIRECT-minVar'
for polIdx = 1:length(selPolicyList)

    selPolicy = selPolicyList{polIdx};
    fprintf('Executing with policy "%s"\n',selPolicy);
    
    if strcmp(selPolicy,'random');  iterAverage = 1;  % Random executed several times and averaged
    else;                           iterAverage = nIter;  % Any other policy only executed once
    end

    % This loops will be executed only once for 
    while(iterAverage<=nIter)
    
        % Initialize Sampling
        if strcmp(selPolicy,'random') || strcmp(selPolicy,'UM') || strcmp(selPolicy,'EI') || strcmp(selPolicy,'PI')
            % Initialize random space sampling
            combo = combvec(elevList_int,azimList_int);  % Truth table
            samp_ids = randperm(length(combo));  % Generate Random sequence (no EI for now)
            samples = combo(:,samp_ids(1:n_samples));  % Extract samples Elevations and Azimuths
            elev_sample_random = samples(1,:);
            azim_sample_random = samples(2,:);
            elev_sample1 = elev_sample_random(1:minSamples);  % Store the sampled elevation angles
            azim_sample1 = azim_sample_random(1:minSamples);  % Store the sampled azimuth angles
        elseif strcmp(selPolicy,'DIRECT-rand') || strcmp(selPolicy,'DIRECT-minVar')
            % Initialize DIRECT
            selPoss = CBG_DIRECT([0 90],[0 180],minSamples);
            % Plot Initial DIRECT configuration
            if plotFlag;    CBG_DIRECTplot(selPoss,2);   end
            % First, select the samples to try (special case, more than one)
            elev_sample1 = selPoss(1,1,:);  elev_sample1 = elev_sample1(:).';
            azim_sample1 = selPoss(1,2,:);  azim_sample1 = azim_sample1(:).';
            % Second, initialize massively DIRECT (special case, more than one)
            for selTrial = 1:size(selPoss,3)
                selCenter = selPoss(1,:,selTrial);
                selLimElev = selPoss(2,:,selTrial);
                selLimAzim = selPoss(3,:,selTrial);
                % Second, break down the new area into new sub-areas using DIRECT
                addSelPoss = CBG_DIRECT(selLimElev,selLimAzim,newAddPoss);
                % Append new sub-area to global area
                selPoss = cat(3, selPoss, addSelPoss);
            end
            % Remove explored possibility from selection set
            selPoss(:,:,1:minSamples) = [];
        end

        count = 1;  % To index results into storing
        for nSamples = nSamplesList

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
            [pred_gain, max_idx] = max(Zhat(:));  % Predicted gain
            [idx_pred_maxGain_X,idx_pred_maxGain_Y] = ind2sub(size(Zhat),max_idx);
            pred_elevation = X0(idx_pred_maxGain_X,idx_pred_maxGain_Y);  % Selected Elevation
            pred_azimuth = Y0(idx_pred_maxGain_X,idx_pred_maxGain_Y);  % Selected Azimuth
            pred_gain_real = Z0(idx_pred_maxGain_X,idx_pred_maxGain_Y);  % Real achieved gain

            if plotFlag

                % Plot DIRECT map
                if strcmp(selPolicy,'DIRECT-rand') || strcmp(selPolicy,'DIRECT-minVar')
                    CBG_DIRECTplot(selPoss,2);
                end

                % Plot exhaustive, predictions and uncertainty
%                 figure(1); 
% 
%                 ax = subplot(1,4,1); cla reset; hold on;
                ax = figure;  hold on;
%                 brewermap('*Spectral');  colormap(ax,'brewermap');  colorbar('eastoutside');
%                 colormap(flipud(gray(15)))
%                 colormap(flipud(jet(15)))
                colormap(jet(15))
                colorbar('eastoutside');
                colorbar('hide')
                imagesc(X0(1,:),Y0(:,1),Z0); axis xy; axis tight; % Generate Exhaustive image
                p1 = plot(elev_sample1,azim_sample1,'o','color','b','MarkerSize',5,'MarkerFaceColor','b');  % Plot trials so far
                p2 = plot(exhv_elevation,exhv_azimuth,'sr','MarkerSize',7.5,'MarkerFaceColor','r','MarkerEdgeColor','k');  % Plot trials so far
                p3 = plot(pred_elevation,pred_azimuth,'sg','MarkerSize',7.5,'MarkerFaceColor','g','MarkerEdgeColor','k');  % Plot trials so far
                mystr = strcat('trial',{' '},num2str(nSamples));
                title(mystr);
                if strcmp(selPolicy,'DIRECT-rand') || strcmp(selPolicy,'DIRECT-minVar')
                    p4 = scatter(selPoss(1,1,:),selPoss(1,2,:),15,'Marker','+','MarkerFaceColor','k','MarkerEdgeColor','k');
                    for t = 1:size(selPoss,3)
                    hm = plot([selPoss(2,1,t) selPoss(2,1,t)],[selPoss(3,1,t) selPoss(3,2,t)],'lineWidth',0.5,'color','k','lineStyle','-');  hm.Color(4)=0.7;
                    hm = plot([selPoss(2,1,t) selPoss(2,2,t)],[selPoss(3,1,t) selPoss(3,1,t)],'lineWidth',0.5,'color','k','lineStyle','-');  hm.Color(4)=0.7;
                    end
                    legend([p1 p2 p3 p4],{'Trials','Optimum','Selected','Candidates'},'Location','NorthWest')
                else
                    legend([p1 p2 p3],{'Trials','Optimum','Selected'},'Location','NorthWest')
                end
                yticks([0 30 60 90 120 150 180]);
                xticks([0 15 30 45 60 75 90]);
                pos = get(gcf, 'Position');
                set(gcf,'position',[pos(1),pos(2),323,244]);
                hold off

%                 ax = subplot(1,4,2); cla reset;
                ax = figure;
                imagesc(X0(1,:),Y0(:,1),Zhat); axis xy; axis tight;
                brewermap('*Spectral');
                colormap(ax,'brewermap');
                colorbar('eastoutside');
                mystr = strcat('#',num2str(nSamples),{' '},'kriging predictions');
                title(mystr);

%                 ax = subplot(1,4,3); cla reset;
                ax = figure; hold on; 
%                 brewermap('*Spectral');  colormap(ax,'brewermap');  colorbar('eastoutside');
                colormap(jet(15))
                colorbar('eastoutside');
                colorbar('hide')
                contourf(X0,Y0,Zvar);
                mystr = strcat('#',num2str(nSamples),{' '},'kriging variance');
                title(mystr);
                yticks([0 30 60 90 120 150 180]);
                xticks([0 15 30 45 60 75 90]);
                pos = get(gcf, 'Position');
                set(gcf,'position',[pos(1),pos(2),323,244]);
                p1 = plot(elev_sample1,azim_sample1,'o','color','g','MarkerSize',5,'MarkerFaceColor','g');  % Plot trials so far
                if strcmp(selPolicy,'DIRECT-rand') || strcmp(selPolicy,'DIRECT-minVar')
                    scatter(selPoss(1,1,:),selPoss(1,2,:),15,'Marker','+','MarkerFaceColor','k','MarkerEdgeColor','k');
                    for t = 1:size(selPoss,3)
                    hm = plot([selPoss(2,1,t) selPoss(2,1,t)],[selPoss(3,1,t) selPoss(3,2,t)],'lineWidth',0.3,'color','k','lineStyle','-');  hm.Color(4)=1;
                    hm = plot([selPoss(2,1,t) selPoss(2,2,t)],[selPoss(3,1,t) selPoss(3,1,t)],'lineWidth',0.3,'color','k','lineStyle','-');  hm.Color(4)=1;
                    end
                end
                hold off

%                 ax = subplot(1,4,4); cla reset;
                ax = figure; 
                brewermap('*RdBu');  colormap(ax,'brewermap');  colorbar('eastoutside');
                plot(abs(Zvar(:)));
                mystr = strcat('#',num2str(nSamples),{' '},'Overall variance');
                title(mystr);
            end

            % Select next trial upong selection policy
            if strcmp(selPolicy,'random')
                elev_sample1 = elev_sample_random(1:nSamples);
                azim_sample1 = azim_sample_random(1:nSamples);
                % Store results - variance
                avVar_random1(count,iterAverage) = mean(abs(Zvar(:)));
                maxVar_random1(count,iterAverage) = max(abs(Zvar(:)));
                minVar_random1(count,iterAverage) = min(abs(Zvar(:)));
                % Store results - predictions
                ml_random_pred_azim(count,iterAverage) = pred_azimuth;
                ml_random_pred_elev(count,iterAverage) = pred_elevation;
                ml_random_pred_gain1(count,iterAverage) = pred_gain;
                ml_random_real_gain1(count,iterAverage) = pred_gain_real;
            elseif strcmp(selPolicy,'EI')
                % To-do
            elseif strcmp(selPolicy,'PI')
                % First, select sample with lowest variance
                [~,trialsSort] = sort(Zhat(:),'descend') ;
                selTrials = trialsSort(1:new_additions).';
                for selTrial = selTrials
                    % Something here
                    [idx_exhv_maxGain_X,idx_exhv_maxGain_Y] = ind2sub(size(Zhat),selTrial);
                    sel_elevation = X0(idx_exhv_maxGain_X,idx_exhv_maxGain_Y);
                    sel_azimuth = Y0(idx_exhv_maxGain_X,idx_exhv_maxGain_Y);
                    % Select actual azimuth and elevation to evaluate with kriging
                    elev_sample1 = [elev_sample1 sel_elevation];  %#ok<AGROW>
                    azim_sample1 = [azim_sample1 sel_azimuth];  %#ok<AGROW>
                end
                % Store results - variance
                avVar_random1(count,iterAverage) = mean(abs(Zvar(:)));
                maxVar_random1(count,iterAverage) = max(abs(Zvar(:)));
                minVar_random1(count,iterAverage) = min(abs(Zvar(:)));
                % Store results - predictions
                ml_random_pred_azim(count,iterAverage) = pred_azimuth;
                ml_random_pred_elev(count,iterAverage) = pred_elevation;
                ml_random_pred_gain1(count,iterAverage) = pred_gain;
                ml_random_real_gain1(count,iterAverage) = pred_gain_real;
            elseif strcmp(selPolicy,'UM')
                % First, select sample with lowest variance
                [~,trialsSort] = sort(Zvar(:),'descend') ;
                selTrials = trialsSort(1:new_additions).';
                for selTrial = selTrials
                    % Something here
                    [idx_exhv_maxGain_X,idx_exhv_maxGain_Y] = ind2sub(size(Zvar),selTrial);
                    sel_elevation = X0(idx_exhv_maxGain_X,idx_exhv_maxGain_Y);
                    sel_azimuth = Y0(idx_exhv_maxGain_X,idx_exhv_maxGain_Y);
                    % Select actual azimuth and elevation to evaluate with kriging
                    elev_sample1 = [elev_sample1 sel_elevation];  %#ok<AGROW>
                    azim_sample1 = [azim_sample1 sel_azimuth];  %#ok<AGROW>
                end
                % Store results - variance
                avVar_UM(count) = mean(abs(Zvar(:)));
                maxVar_UM(count) = max(abs(Zvar(:)));
                minVar_UM(count) = min(abs(Zvar(:)));
                % Store results - predictions
                ml_UM_pred_azim(count) = pred_azimuth;
                ml_UM_pred_elev(count) = pred_elevation;
                ml_UM_pred_gain(count) = pred_gain;
                ml_UM_real_gain(count) = pred_gain_real;
                ml_UM_pred_gap_gain(count) = exhv_gain - pred_gain;
                ml_UM_real_gap_gain(count) = exhv_gain - pred_gain_real;
            elseif strcmp(selPolicy,'DIRECT-rand')
                % First, select sample following a random selection
                selTrials = randperm(size(selPoss,3));
                selTrials = selTrials(1:new_additions);
                for selTrial = selTrials
                    % Retrieve results
                    selCenter = selPoss(1,:,selTrial);
                    selLimElev = selPoss(2,:,selTrial);
                    selLimAzim = selPoss(3,:,selTrial);
                    % Second, break down the new area into new sub-areas using DIRECT
                    addSelPoss = CBG_DIRECT(selLimElev,selLimAzim,newAddPoss);
                    % Append new sub-area to global area
                    selPoss = cat(3, selPoss, addSelPoss);
                    % Select actual azimuth and elevation to evaluate with kriging
                    elev_sample1 = [elev_sample1 selCenter(1)];  %#ok<AGROW>
                    azim_sample1 = [azim_sample1 selCenter(2)];  %#ok<AGROW>
                end
                % Remove explored possibility from selection set
                selPoss(:,:,selTrials) = [];
                % Store results - variance
                avVar_DIRECTrandom(count) = mean(abs(Zvar(:)));
                maxVar_DIRECTrandom(count) = max(abs(Zvar(:)));
                minVar_DIRECTrandom(count) = min(abs(Zvar(:)));
                % Store results - predictions
                ml_DIRECTrandom_pred_azim(count) = pred_azimuth;
                ml_DIRECTrandom_pred_elev(count) = pred_elevation;
                ml_DIRECTrandom_pred_gain(count) = pred_gain;
                ml_DIRECTrandom_real_gain(count) = pred_gain_real;
                ml_DIRECTrandom_pred_gap_gain(count) = exhv_gain - ml_DIRECTrandom_pred_gain(count);
                ml_DIRECTrandom_real_gap_gain(count) = exhv_gain - ml_DIRECTrandom_real_gain(count);
            elseif strcmp(selPolicy,'DIRECT-minVar')
                % First, select sample with lowest variance
                myElevs = selPoss(1,1,:);   myElevs = myElevs(:);
                myAzims = selPoss(1,2,:);   myAzims = myAzims(:);
                myX0 = X0(:);  myY0 = Y0(:);  myZvar = Zvar(:);
                myX01 = repmat(myX0,1,length(myAzims)).';
                myY01 = repmat(myY0,1,length(myElevs)).';
                myIndices = zeros(1,length(myAzims));
                for temp = 1:length(myAzims)
                    test_azim = find(myY0==ceil(myAzims(temp)));  %Swap them, needed
                    test_elev = find(myX0==ceil(myElevs(temp)));  %Swap them, needed
                    myIndices(temp) = intersect(test_azim,test_elev);
                end
                [B,trialsSort] = sort(myZvar(myIndices),'descend') ;
                selTrials = trialsSort(1:new_additions).';
                for selTrial = selTrials
                    % Retrieve results
                    selCenter = selPoss(1,:,selTrial);
                    selLimElev = selPoss(2,:,selTrial);
                    selLimAzim = selPoss(3,:,selTrial);
                    % Second, break down the new area into new sub-areas using DIRECT
                    addSelPoss = CBG_DIRECT(selLimElev,selLimAzim,newAddPoss);
                    % Append new sub-area to global area
                    selPoss = cat(3, selPoss, addSelPoss);
                    % Select actual azimuth and elevation to evaluate with kriging
                    elev_sample1 = [elev_sample1 selCenter(1)];  %#ok<AGROW>
                    azim_sample1 = [azim_sample1 selCenter(2)];  %#ok<AGROW>
                end
                % Remove explored possibility from selection set
                selPoss(:,:,selTrials) = [];
                % Store results - variance
                avVar_DIRECTminVar(count) = mean(abs(Zvar(:)));
                maxVar_DIRECTminVar(count) = max(abs(Zvar(:)));
                minVar_DIRECTminVar(count) = min(abs(Zvar(:)));
                % Store results - predictions
                ml_DIRECTminVar_pred_azim(count) = pred_azimuth;
                ml_DIRECTminVar_pred_elev(count) = pred_elevation;
                ml_DIRECTminVar_pred_gain(count) = pred_gain;
                ml_DIRECTminVar_real_gain(count) = pred_gain_real;
                ml_DIRECTminVar_pred_gap_gain(count) = exhv_gain - pred_gain;
                ml_DIRECTminVar_real_gap_gain(count) = exhv_gain - pred_gain_real;
            end

            % Update counter for results
            count = count + 1;
        end

        % Average Random results and store in global variable
        if strcmp(selPolicy,'random')
            % Store results - variance
            avVar_random = mean(avVar_random1,2);
            maxVar_random = mean(maxVar_random1,2);
            minVar_random = mean(minVar_random1,2);
            % Store results - gain
            ml_random_pred_gain = mean(ml_random_pred_gain1,2);
            ml_random_real_gain = mean(ml_random_real_gain1,2);
            ml_random_pred_gap_gain = exhv_gain - ml_random_pred_gain;
            ml_random_real_gap_gain = exhv_gain - ml_random_real_gain;
        end
    
        % Increase averaging index
        iterAverage = iterAverage + 1;
    end
end

figIdx = 50;
%% AVERAGE VARIANCE
% figure(figIdx); figIdx = figIdx + 1; hold on;
% plot(nSamplesList,avVar_random,'lineWidth',2);
% plot(nSamplesList,avVar_DIRECTrandom,'lineWidth',2);
% plot(nSamplesList,avVar_DIRECTminVar,'lineWidth',2);
% grid minor;
% lg = legend('Random - Random','Direct - Random sel.','Direct - minVar sel.','Optimum');
% set(lg,'FontSize',10);
% title('analysis on the AVERAGE Variance achieved','FontSize',12);
% xlabel('Number of trials','FontSize',12);
% ylabel('Variance','FontSize',12);
% hold off;
%% MAXIMUM VARIANCE
% figure(figIdx); figIdx = figIdx + 1; hold on;
% plot(nSamplesList,maxVar_random,'lineWidth',2);
% plot(nSamplesList,maxVar_DIRECTrandom,'lineWidth',2);
% plot(nSamplesList,maxVar_DIRECTminVar,'lineWidth',2);
% grid minor;
% lg = legend('Random - Random','Direct - Random sel.','Direct - minVar sel.','Optimum');
% set(lg,'FontSize',10);
% title('analysis on the MAX Variance achieved','FontSize',12);
% xlabel('Number of trials','FontSize',12);
% ylabel('Variance','FontSize',12);
%% MINIMUM VARIANCE
% figure(figIdx); figIdx = figIdx + 1; hold on;
% plot(nSamplesList,minVar_random,'lineWidth',2);
% plot(nSamplesList,minVar_DIRECTrandom,'lineWidth',2);
% plot(nSamplesList,minVar_DIRECTminVar,'lineWidth',2);
% grid minor;
% lg = legend('Random - Random','Direct - Random sel.','Direct - minVar sel.','Optimum');
% set(lg,'FontSize',10);
% title('analysis on the MIN Variance achieved','FontSize',12);
% xlabel('Number of trials','FontSize',12);
% ylabel('Variance','FontSize',12);
% hold off
%% PREDICTED GAIN
figure(figIdx); figIdx = figIdx + 1; hold on;
plot(nSamplesList,ml_random_pred_gain,'lineWidth',2);
plot(nSamplesList,ml_UM_pred_gain,'lineWidth',2);
plot(nSamplesList,ml_DIRECTrandom_pred_gain,'lineWidth',2);
plot(nSamplesList,ml_DIRECTminVar_pred_gain,'lineWidth',2);
plot(nSamplesList,exhv_gain_long,'lineWidth',2);
grid minor;
lg = legend('Random - Random','Uncertainty min.','Direct - Random sel.','Direct - minVar sel.','Optimum');
set(lg,'FontSize',8);
title('PREDICTED MAX GAIN achieved','FontSize',12);
xlabel('Number of trials','FontSize',12);
ylabel('Variance','FontSize',12);
hold off;
%% PREDICTED GAIN (BAR)
% groupLabels = {nSamplesList};  % labels to use on tick marks for groups
% stackData = [ml_random_pred_gain ml_UM_pred_gain.' ml_DIRECTrandom_pred_gain.' ml_DIRECTminVar_pred_gain.' exhv_gain_long.'];
% plotBarStackGroups(stackData, groupLabels,figIdx);  figIdx = figIdx + 1; hold on; 
% lg = legend('Random - Random','Uncertainty min.','Direct - Random sel.','Direct - minVar sel.','Optimum');
% set(lg,'FontSize',8);
% title('analysis on the MAX GAIN achieved','FontSize',12);
% xlabel('Number of trials','FontSize',12);
% ylabel('channel gain (linear)','FontSize',12);
% grid minor;
%% PREDICTED AZIMUTH
% figure(figIdx); figIdx = figIdx + 1; hold on;
% plot(nSamplesList,ml_random_pred_azim,'lineWidth',2);
% plot(nSamplesList,ml_UM_pred_azim,'lineWidth',2);
% plot(nSamplesList,ml_DIRECTrandom_pred_azim,'lineWidth',2);
% plot(nSamplesList,ml_DIRECTminVar_pred_azim,'lineWidth',2);
% plot(nSamplesList,exhv_azimuth,'lineWidth',2);
% grid minor;
% lg = legend('Random - Random','Uncertainty min.','Direct - Random sel.','Direct - minVar sel.','Optimum');
% set(lg,'FontSize',8);
% title('analysis on the AZIMUTH prediction','FontSize',12);
% xlabel('Number of trials','FontSize',12);
% ylabel('Variance','FontSize',12);
% hold off;
%% PREDICTED ELEVATION
% figure(figIdx); figIdx = figIdx + 1; hold on; 
% plot(nSamplesList,ml_random_pred_elev,'lineWidth',2);
% plot(nSamplesList,ml_UM_pred_elev,'lineWidth',2);
% plot(nSamplesList,ml_DIRECTrandom_pred_elev,'lineWidth',2);
% plot(nSamplesList,ml_DIRECTminVar_pred_elev,'lineWidth',2);
% plot(nSamplesList,exhv_elevation,'lineWidth',2);
% grid minor;
% lg = legend('Random - Random','Uncertainty min.','Direct - Random sel.','Direct - minVar sel.','Optimum');
% set(lg,'FontSize',8);
% title('analysis on the AZIMUTH prediction','FontSize',12);
% xlabel('Number of trials','FontSize',12);
% ylabel('Variance','FontSize',12);
% hold off;
%% REAL GAIN ACHIEVED
figure(figIdx); figIdx = figIdx + 1; hold on;
plot(nSamplesList,ml_random_real_gain,'lineWidth',2);
plot(nSamplesList,ml_UM_real_gain,'lineWidth',2);
plot(nSamplesList,ml_DIRECTrandom_real_gain,'lineWidth',2);
plot(nSamplesList,ml_DIRECTminVar_real_gain,'lineWidth',2);
plot(nSamplesList,exhv_gain_long,'lineWidth',2);
grid minor;
lg = legend('Random - Random','Uncertainty min.','Direct - Random sel.','Direct - minVar sel.','Optimum');
set(lg,'FontSize',8);
title('REAL MAX GAIN achieved','FontSize',12);
xlabel('Number of trials','FontSize',12);
ylabel('Variance','FontSize',12);
hold off;
%% REAL GAIN (BAR)
% groupLabels = {nSamplesList};  % labels to use on tick marks for groups
% stackData = [ml_random_real_gain ml_DIRECTrandom_real_gain.' ml_DIRECTminVar_real_gain.' exhv_gain_long.'];
% plotBarStackGroups(stackData, groupLabels,figIdx);  figIdx = figIdx + 1; hold on; 
% lg = legend('RANDOM-samp and RANDOM-sel','DIRECT-samp and RANDOM-sel','DIRECT-samp and MINVAR-sel','Optimum Exhaustive');
% set(lg,'FontSize',12);
% grid minor;
%% REAL GAP TO OPTIMUM
figure(figIdx); figIdx = figIdx + 1; hold on;
plot(nSamplesList,ml_random_real_gap_gain,'lineWidth',2);
plot(nSamplesList,ml_UM_real_gap_gain,'lineWidth',2);
plot(nSamplesList,ml_DIRECTrandom_real_gap_gain,'lineWidth',2);
plot(nSamplesList,ml_DIRECTminVar_real_gap_gain,'lineWidth',2);
grid minor;
lg = legend('Random - Random','Uncertainty min.','Direct - Random sel.','Direct - minVar sel.');
set(lg,'FontSize',8);
title('GAP TO REAL MAX GAIN','FontSize',12);
xlabel('Number of trials','FontSize',12);
ylabel('Variance','FontSize',12);
hold off;