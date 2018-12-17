function [opt_krig,opt_exhv,Zhat,Zvar] = CBG_kriging_old(Z,elevList,azimList,n_samples,plotFlag)

nTxAntennas = size(Z,3);
[X,Y] = meshgrid(elevList,azimList);
plotIdx = 1;

% Output variables
opt_krig = zeros(3,nTxAntennas);  % Kriging -  gain, elevation and azimuth (deg)
opt_exhv = zeros(3,nTxAntennas);  % Exhaustive - gain, elevation and azimuth (deg)

for id = 1:nTxAntennas
    % New interpolated input (assume this is the real exhaustive)
    elevList_int = (max(elevList):-5:min(elevList));
    azimList_int = (max(azimList):-5:min(azimList));
    [X0,Y0] = meshgrid(elevList_int,azimList_int);
    Z0(:,:,id) = interp2(X,Y,Z(:,:,id),X0,Y0,'spline');  %#ok

    % Sample the space
    combo = combvec(elevList_int,azimList_int);  % Truth table
    samp_ids = randperm(length(combo));  % Generate Random sequence (no EI for now)
    samples = combo(:,samp_ids(1:n_samples));  % Extract samples Elevations and Azimuths
    elev_sample = samples(1,:);
    azim_sample = samples(2,:);
    [X_sample,Y_sample] = meshgrid(elev_sample,azim_sample);
    Z_sample = interp2(X0,Y0,Z0(:,:,id),X_sample,Y_sample);

    % Calculate the sample variogram
    a = reshape(X_sample,[size(X_sample,1)*size(X_sample,2),1]);
    b = reshape(Y_sample,[size(Y_sample,1)*size(Y_sample,2),1]);
    c = reshape(Z_sample,[size(Z_sample,1)*size(Z_sample,2),1]);
    v = variogram([a b],c,'plotit',false,'maxdist',100);

    % Compute Gaussian fit following Least Squares
    [~,~,~,vstruct] = variogramfit(v.distance,v.val,[],[],[],'model','gaussian','plotit',false);

    % Now use the sampled locations in a kriging
    [Zhat(:,:,id),Zvar(:,:,id)] = kriging(vstruct,X_sample,Y_sample,Z_sample,X0,Y0);  %#ok
    
    % Find Maximum based on prediction and uncertainty
    T = Zhat(:,:,id);
    [opt_krig(1,id),B] = max(T(:));
    [idx_azim,idx_elev] = ind2sub(size(T),B);
    opt_krig(2,id) = elevList_int(idx_elev);
    opt_krig(3,id) = azimList_int(idx_azim);
    T = Z(:,:,id);
    [opt_exhv(1,id),B] = max(T(:));
    [idx_azim,idx_elev] = ind2sub(size(T),B);
    opt_exhv(2,id) = elevList_int(idx_elev);
    opt_exhv(3,id) = azimList_int(idx_azim);

    % Plot results
    if plotFlag
        % Generate exhaustive image
        subplot(nTxAntennas,4,plotIdx); hold on;  plotIdx = plotIdx + 1;
        imagesc(X0(1,:),Y0(:,1),Z0(:,:,id)); axis image; axis xy
        plot(elev_sample,azim_sample,'+r','MarkerSize',5);
        title('Real Map and sampling');
        
        subplot(nTxAntennas,4,plotIdx); hold on;  plotIdx = plotIdx + 1;
        variogramfit(v.distance,v.val,[],[],[],'model','gaussian','plotit',true);
        title('Prior fitting: Spatial correlation');
        
        subplot(nTxAntennas,4,plotIdx); hold on;  plotIdx = plotIdx + 1;
        imagesc(X0(1,:),Y0(:,1),Zhat(:,:,id)); axis image; axis xy
        title('kriging predictions');
    
        subplot(nTxAntennas,4,plotIdx); hold on;  plotIdx = plotIdx + 1;
        contour(X0,Y0,Zvar(:,:,id)); axis image
        title('kriging variance')
    end
end


% EOF