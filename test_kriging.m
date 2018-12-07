% clear all; clear classes;  %#ok
close all; clc;
addpath('kriging/');  % Include Kriging folder (variogram and kriging)
addpath('BrewerMap/');  % Include additional colors: 'BrBG'|'PRGn'|'PiYG'|'PuOr'|'RdBu'|'RdGy'|'RdYlBu'|'RdYlGn'|'Spectral'

% load('results/res_ROT_outdoors.mat','azimList','elevList','chTot','nTxAntennas','appliedElev','appliedAzym');
% load('results-19-Nov-2018_16:17:56.mat','azimList','elevList','chTot','nTxAntennas','appliedElev','appliedAzym');
% load('results-19-Nov-2018_18:52:28.mat','azimList','elevList','chTot','nTxAntennas','appliedElev','appliedAzym');
% load('results-19-Nov-2018_18:52:28_2.mat','azimList','elevList','chTot','nTxAntennas','appliedElev','appliedAzym');
% load('results-25-Nov-2018_14:49:13.mat','azimList','elevList','chTot','nTxAntennas','appliedElev','appliedAzym');
% load('results-25-Nov-2018_14:55:08.mat','azimList','elevList','chTot','nTxAntennas','appliedElev','appliedAzym');
% load('results-25-Nov-2018_20:15:13.mat','azimList','elevList','chTot','nTxAntennas','appliedElev','appliedAzym');

% load('results-27-Nov-2018_01:32:43.mat','azimList','elevList','chTot','nTxAntennas','appliedElev','appliedAzym');

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
n_samples = 20;  % Maximum number of samples to run Kriging with
minSamples = 3;  % Minimum Kriging initialization space
new_additions = 1;
id = 2; % For instance, could be 1,2,3,4 (antenna id)

% New interpolated input (assume this is the real exhaustive)
elevList_int = (max(elevList):-1:min(elevList));
azimList_int = (max(azimList):-1:min(azimList));
[X0,Y0] = meshgrid(elevList_int,azimList_int);
Z0 = interp2(X,Y,Z(:,:,id),X0,Y0,'spline');

% Generate Exhaustive image

%     % Combo 1    
%     elev_sample = [20 20 20 60 60 60 80 80 80];
%     azim_sample = [40 90 140 40 90 140 40 90 140];
%     % Combo 2
%     elev_sample = [20 20 80 80];
%     azim_sample = [40 60 40 60];
%     % Combo 3
%     elev_sample = [20 20 40 40];
%     azim_sample = [40 140 40 140];

% Sample the space
combo = combvec(elevList_int,azimList_int);  % Truth table
samp_ids = randperm(length(combo));  % Generate Random sequence (no EI for now)
samples = combo(:,samp_ids(1:n_samples));  % Extract samples Elevations and Azimuths
elev_sample = samples(1,:);
azim_sample = samples(2,:);

count = 1;
myvar = zeros(size(minSamples:new_additions:n_samples));
for iter = minSamples:new_additions:n_samples
    elev_sample1 = elev_sample(1:iter);
    azim_sample1 = azim_sample(1:iter);

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
    
    % plotting
	ax = subplot(1,4,1); hold on
    brewermap('*RdBu');
    colormap(ax,'brewermap');
    colorbar('eastoutside');
    imagesc(X0(1,:),Y0(:,1),Z0); axis image; axis xy
    plot(elev_sample1,azim_sample1,'+r','MarkerSize',5);
    mystr = strcat('#',num2str(iter),{' '},'Real+sampling');
    title(mystr);
    hold off

    ax = subplot(1,4,2);
    imagesc(X0(1,:),Y0(:,1),Zhat); axis image; axis xy;
    brewermap('*RdBu');
    colormap(ax,'brewermap');
    colorbar('eastoutside');
    mystr = strcat('#',num2str(iter),{' '},'kriging predictions');
    title(mystr);

    ax = subplot(1,4,3);
    contourf(X0,Y0,Zvar);
    mystr = strcat('#',num2str(iter),{' '},'kriging variance');
    title(mystr);
    brewermap('*RdBu'); % 'BrBG'|'PRGn'|'PiYG'|'PuOr'|'RdBu'|'RdGy'|'RdYlBu'|'RdYlGn'|'Spectral'
    colormap(ax,'brewermap');
    colorbar('eastoutside');
    
    ax = subplot(1,4,4);
    A = reshape(Zvar,size(Zvar,1)*size(Zvar,2),1);
    plot(abs(A));
    brewermap('*RdBu');
    colormap(ax,'brewermap');
    colorbar('eastoutside');
    mystr = strcat('#',num2str(iter),{' '},'Overall variance');
    title(mystr);
    
    myvar(count) = max(abs(A));
    count = count + 1;
end

figure;
plot(minSamples:new_additions:n_samples,myvar);