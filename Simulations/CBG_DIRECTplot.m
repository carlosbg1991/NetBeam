function ax = CBG_DIRECTplot(selPoss,figIdx)
% CBG_DIRECTplot - Plots the Rectangle area to be explored
%
% Syntax:  ax = CBG_DIRECTplot(selPoss,figIdx)
%
% Inputs:
%    selPoss - Matrix containing information of centers and limits for each
%              point
%    figIdx - Figure index for the plot
%
% Outputs:
%    ax - Figure handle of the generated plot
%
% Example: 
%    To-do
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: To-do

%------------- BEGIN CODE --------------
ax = figure(figIdx); cla reset; hold on;
scatter(selPoss(1,1,:),selPoss(1,2,:),15,'Marker','s','MarkerFaceColor',[0 .7 .7],'MarkerEdgeColor',[0 .5 .5]);
for t = 1:size(selPoss,3)
    plot([selPoss(2,1,t) selPoss(2,1,t)],[selPoss(3,1,t) selPoss(3,2,t)],'lineWidth',1.2,'color','k','lineStyle','-');
    plot([selPoss(2,1,t) selPoss(2,2,t)],[selPoss(3,1,t) selPoss(3,1,t)],'lineWidth',1.2,'color','k','lineStyle','-');
end
hold off;


%EOF