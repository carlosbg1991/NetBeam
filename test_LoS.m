% clear all; clear classes;  %#ok
close all; clc;
addpath('BrewerMap/');  % Include additional colors: 'BrBG'|'PRGn'|'PiYG'|'PuOr'|'RdBu'|'RdGy'|'RdYlBu'|'RdYlGn'|'Spectral'
addpath('export_fig/');  % export figure in eps

%% PARAMETERS
antIDList  = (1:4);     % Antenna ID, could be 1,2,3,4
expIDList  = (1:5);     % Experiment ID, could be 1,2,3,4,5

%% PARSE Data if not done before
if ~exist('outdoor','var') || ~exist('indoor','var')
    CBG_parse_experiments;  % Parse data into structs 'indoor', 'outdoor'
end
% indoor = 2;  % useful to change variable name
elevList = indoor.elevList;
azimList = indoor.azimList;

for expID = expIDList
    for antID = antIDList
        [X,Y] = meshgrid(elevList,azimList);  % Orientation space (exhaustive)
        Z = indoor.gainTot(:,:,antID,expID);  % Channel gain (exhaustive)

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
        los_elevation = indoor.elevLoS(expID,antID);
        los_azimuth = indoor.azimLoS(expID,antID);
        
        ax = figure;  hold on;
        
%         axes1 = axes('Parent',ax);  hold(axes1,'on');
%         p1 = surface(X0,Y0,Z0,'lineStyle','none');
        imagesc(X0(1,:),Y0(:,1),Z0); axis xy; axis tight; % Generate Exhaustive image
        p2 = plot(exhv_elevation,exhv_azimuth,'sk','MarkerSize',10,'MarkerFaceColor','k','MarkerEdgeColor','k');  % Plot trials so far
        p3 = plot(los_elevation,los_azimuth,'sg','MarkerSize',10,'MarkerFaceColor','g','MarkerEdgeColor','k');  % Plot trials so far
        lg = legend([p2 p3],'Max. gain','LoS');
        h = get(gca,'Children');
        set(gca,'Children',[h(1) h(2) h(3)]);
        colormap(jet(15))
        colorbar('eastoutside');
        colorbar('hide')
        hold off
        axis tight; 
        
%         h = get(gca,'Children');
%         set(gca,'Children',[h(3) h(2) h(1)])
        
%         % Create colorbar
%         colorbar(axes1);
%         axis(axes1,'tight');
%         % Set the remaining axes properties
%         set(axes1,'CLim',[0.01619 0.61566],'Colormap',...
%             [0.650980412960052 0.650980412960052 0.650980412960052;0.650980412960052 0.650980412960052 0.650980412960052;0.488235294818878 0.688235282897949 0.738235294818878;0.325490206480026 0.725490212440491 0.825490236282349;0.162745103240013 0.762745141983032 0.912745118141174;0 0.800000011920929 1;0.25 0.850000023841858 0.75;0.5 0.899999976158142 0.5;0.75 0.949999988079071 0.25;1 1 0;1 0.75 0;1 0.5 0;1 0.25 0;1 0 0;0.75 0 0]);
%         % Create colorbar
%         colorbar(axes1);

        % Try this size for including it in prel. results in the paper
        pos = get(gcf, 'Position');
        set(gcf,'position',[pos(1),pos(2),250,250]);
        yticks([0 30 60 90 120 150 180]);
        xticks([0 15 30 45 60 75 90]);
        
    end
end