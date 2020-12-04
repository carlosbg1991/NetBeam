% clear all; clear classes;  %#ok
close all; clc;
addpath('../BrewerMap/');  % Include additional colors: 'BrBG'|'PRGn'|'PiYG'|'PuOr'|'RdBu'|'RdGy'|'RdYlBu'|'RdYlGn'|'Spectral'
addpath('kriging/');  % Include variogram and variogramfit
addpath('data/');  % where results from experiments using real radios are

%% PARAMETERS
antIDList  = (1:3);     % Antenna ID, could be 1,2,3,4
expIDList  = (1:3);     % Experiment ID, could be 1,2,3,4,5
colOrder = colororder;

%% PARSE Data if not done before
if ~exist('RESULTS','var')
    load('RESULTS','indoor','outdoor','paramList');
elseif ~exist('outdoor','var') || ~exist('indoor','var')
    CBG_parse_experiments;  % parse experimental DATA
end
% indoor = 1;
elevList = indoor.elevList;
azimList = indoor.azimList;
[X,Y] = meshgrid(elevList,azimList);  % Orientation space (exhaustive)

idx = 1;
hp = [];
for expID = expIDList
    for antID = antIDList
        % get color ID
        colId = mod(((expID-1)*length(antIDList)) + antID,length(colOrder)) + 1;
        % Channel gain (exhaustive)
        Z = indoor.gainTot(:,:,antID,expID);
        % calculate the sample variogram
        a = reshape(X,[size(X,1)*size(X,2),1]);
        b = reshape(Y,[size(Y,1)*size(Y,2),1]);
        c = reshape(Z,[size(Z,1)*size(Z,2),1]);
        v = variogram([a b],c,'maxdist',150,'plotit',false);
        [~,~,~,vstruct] = variogramfit(v.distance,v.val,[],[],[],'model','gaussian','plotit',false);
        figure(1);  hold on
        p = plot(vstruct.h,vstruct.gammahat,'color',colOrder(colId,:),'lineWidth',2);
        hp1 = plot(v.distance,v.val,'.','lineStyle','none','lineWidth',2,'MarkerSize',15,'color',colOrder(colId,:));
        hp = [hp hp1];
        myString(idx) = strcat('Exp. #',num2str(expID),{' '},'Antenna #',num2str(antID));
        idx = idx + 1;
    end
end
myString{idx} = 'Gaussian fit';
legend([hp p],myString);
ylabel('\gamma(h)','FontSize',12)
xlabel('Angular distance','FontSize',12)
title('Indoor - Gaussian fit using BLUE','FontSize',12)
grid minor