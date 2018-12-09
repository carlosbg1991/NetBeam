% clear all; clear classes;  %#ok
close all; clc;
addpath('BrewerMap/');  % Include additional colors: 'BrBG'|'PRGn'|'PiYG'|'PuOr'|'RdBu'|'RdGy'|'RdYlBu'|'RdYlGn'|'Spectral'

%% PARAMETERS
antIDList  = (1:4);     % Antenna ID, could be 1,2,3,4
expIDList  = (1:5);     % Experiment ID, could be 1,2,3,4,5

%% PARSE Data if not done before
if ~exist('outdoor','var') || ~exist('indoor','var')
    CBG_parse_experiments;  % Parse data into structs 'indoor', 'outdoor'
end
% outdoor = 2;  % useful to change variable name
elevList = outdoor.elevList;
azimList = outdoor.azimList;

for expID = expIDList
    for antID = antIDList
        [X,Y] = meshgrid(elevList,azimList);  % Orientation space (exhaustive)
        Z = outdoor.gainTot(:,:,antID,expID);  % Channel gain (exhaustive)

        % new interpolated input (assume this is the real exhaustive)
        elevList_int = (max(elevList):-1:min(elevList));
        azimList_int = (max(azimList):-1:min(azimList));
        [X0,Y0] = meshgrid(elevList_int,azimList_int);
        Z0 = interp2(X,Y,Z,X0,Y0,'spline');
        
        % find optimum angles - Exhaustive
        [exhv_gain, max_idx] = max(Z0(:));
        [idx_exhv_maxGain_X,idx_exhv_maxGain_Y] = ind2sub(size(Z0),max_idx);
        exhv_elevation = X0(idx_exhv_maxGain_X,idx_exhv_maxGain_Y);
        exhv_azimuth = Y0(idx_exhv_maxGain_X,idx_exhv_maxGain_Y);
        
        % find forecasted angles
        los_elevation = outdoor.elevLoS(expID,antID);
        los_azimuth = outdoor.azimLoS(expID,antID);
        
        ax = figure; hold on;
        brewermap('*RdYlBu');
        colormap(ax,'brewermap');
        colorbar('eastoutside');
        imagesc(X0(1,:),Y0(:,1),Z0); axis image; axis xy  % Generate Exhaustive image
        p2 = plot(exhv_elevation,exhv_azimuth,'sk','MarkerSize',10,'MarkerFaceColor','k','MarkerEdgeColor',[0 0 205]./255);  % Plot trials so far
        p3 = plot(los_elevation,los_azimuth,'sy','MarkerSize',10,'MarkerFaceColor','y','MarkerEdgeColor','k');  % Plot trials so far
        title('');
        legend([p2 p3],'Max. gain','LoS');
        hold off
        
        % Try this size for including it in prel. results in the paper
        pos = get(gcf, 'Position');
        set(gcf,'position',[pos(1),pos(2),width,height]);
    end
end
