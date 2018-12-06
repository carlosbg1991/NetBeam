clear all; clear classes; close all; clc;  %#ok
addpath('kriging/');  % Include Kriging folder (variogram and kriging)

% load('results/res_ROT_outdoors.mat','azimList','elevList','chTot','nTxAntennas','appliedElev','appliedAzym');
% load('results-19-Nov-2018_16:17:56.mat','azimList','elevList','chTot','nTxAntennas','appliedElev','appliedAzym');
% load('results-19-Nov-2018_18:52:28.mat','azimList','elevList','chTot','nTxAntennas','appliedElev','appliedAzym');
% load('results-19-Nov-2018_18:52:28_2.mat','azimList','elevList','chTot','nTxAntennas','appliedElev','appliedAzym');
% load('results-25-Nov-2018_14:49:13.mat','azimList','elevList','chTot','nTxAntennas','appliedElev','appliedAzym');
load('results-25-Nov-2018_14:55:08.mat','azimList','elevList','chTot','nTxAntennas','appliedElev','appliedAzym');
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

n_samples = 50;
for id = 1:nTxAntennas
    % New interpolated input (assume this is the real exhaustive)
    elevList_int = (max(elevList):-1:min(elevList));
    azimList_int = (max(azimList):-1:min(azimList));
    [X0,Y0] = meshgrid(elevList_int,azimList_int);
    Z0(:,:,id) = interp2(X,Y,Z(:,:,id),X0,Y0,'spline');

    % Generate Exhaustive image
    subplot(nTxAntennas,4,plotIdx); hold on;  plotIdx = plotIdx + 1;
    imagesc(X0(1,:),Y0(:,1),Z0(:,:,id)); axis image; axis xy

    % Sample the space
    indeces_elev = randperm(length(elevList_int));  % indices to sample elevation
    indeces_elev = indeces_elev(1:n_samples);
    indeces_azim = randperm(length(azimList_int));  % indices to sample azimuth
    indeces_azim = indeces_azim(1:n_samples);
    elev_sample = elevList_int(indeces_elev);  % select sample elevation
    azim_sample = azimList_int(indeces_azim);  % select sample azimuth
    [X_sample,Y_sample] = meshgrid(elev_sample,azim_sample);
    Z_sample = interp2(X0,Y0,Z0(:,:,id),X_sample,Y_sample);
    plot(elev_sample,azim_sample,'+r','MarkerSize',5);
    title('Real Map and sampling');

    % calculate the sample variogram
    a = reshape(X_sample,[size(X_sample,1)*size(X_sample,2),1]);
    b = reshape(Y_sample,[size(Y_sample,1)*size(Y_sample,2),1]);
    c = reshape(Z_sample,[size(Z_sample,1)*size(Z_sample,2),1]);
    v = variogram([a b],c,'plotit',false,'maxdist',100);

    % Compute Gaussian fit following Least Squares
    subplot(nTxAntennas,4,plotIdx); hold on;  plotIdx = plotIdx + 1;
    [~,~,~,vstruct] = variogramfit(v.distance,v.val,[],[],[],'model','gaussian');
    title('Prior fitting: Spatial correlation');

    % now use the sampled locations in a kriging
    subplot(nTxAntennas,4,plotIdx); hold on;  plotIdx = plotIdx + 1;
    [Zhat,Zvar] = kriging(vstruct,X_sample,Y_sample,Z_sample,X0,Y0);
    imagesc(X0(1,:),Y0(:,1),Zhat); axis image; axis xy
    title('kriging predictions')

    subplot(nTxAntennas,4,plotIdx); hold on;  plotIdx = plotIdx + 1;
    contour(X0,Y0,Zvar); axis image
    title('kriging variance')
end