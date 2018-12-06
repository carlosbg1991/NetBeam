function [Zhat,Zvar] = CBG_kriging(Z,elevList,azimList,n_samples,plotFlag)

nTxAntennas = size(Z,3);
[X,Y] = meshgrid(elevList,azimList);
plotIdx = 1;

for id = 1:nTxAntennas
    % New interpolated input (assume this is the real exhaustive)
    elevList_int = (max(elevList):-1:min(elevList));
    azimList_int = (max(azimList):-1:min(azimList));
    [X0,Y0] = meshgrid(elevList_int,azimList_int);
    Z0(:,:,id) = interp2(X,Y,Z(:,:,id),X0,Y0,'spline');  %#ok

    % Sample the space
    indeces_elev = randperm(length(elevList_int));  % indices to sample elevation
    indeces_elev = indeces_elev(1:n_samples);
    indeces_azim = randperm(length(azimList_int));  % indices to sample azimuth
    indeces_azim = indeces_azim(1:n_samples);
    elev_sample = elevList_int(indeces_elev);  % select sample elevation
    azim_sample = azimList_int(indeces_azim);  % select sample azimuth
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
    [opt_gain_krig(id),B] = max(T(:));
    [opt_azim_krig(id), opt_elev_krig(id)] = ind2sub(size(T),B);
    T = Z(:,:,id);
    [opt_gain_exhv(id),B] = max(T(:));
    [opt_azim_exhv(id), opt_elev_exhv(id)] = ind2sub(size(T),B);

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