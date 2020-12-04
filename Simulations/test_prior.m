function test_prior(varargin)
% TEST_prior - This scripts validates the Gaussian Kernel used in Kriging
% to model channel correlations in the angular domain. It does so by
% computing the variogram of the empirical data and fitting a Gaussian
% model.
%
% This script assesses Fig. 9 of the publication:
% [1] C. Bocanegra, K. Alemdar, S. Garcia, C. Singhal and K. R. Chowdhury,
%     “NetBeam: Network of Distributed Full-dimension
%     Beamforming SDRs for Multi-user Heterogeneous Traffic,” IEEE Dynamic
%     Spectrum (DySpan), Newark, NJ, 2019
%
% Syntax:  test_prior(environment)
%
% Inputs:
%    environment [optional] - String that takes the values 'indoor' and
%    'outdoor'. It selects the environment in which the experiments were
%    carried out. It defaults to indoor.
%
% Outputs: []
%
%
%------------- BEGIN CODE --------------



if (nargin==1)
    environment = varargin{1};
elseif (nargin==0)
    clear all; clear classes;  %#ok
    close all; clc;
    environment = 'outdoor';  % take indoor for example
else
    error('ERROR: The script only accepts 1 or none inputs\n');
end

fprintf('Selected environment: %s\n',environment);

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

if strcmp(environment, 'outdoor')
    elevList = outdoor.elevList;
    azimList = outdoor.azimList;
    gainTot = outdoor.gainTot;
elseif strcmp(environment, 'indoor')
    elevList = indoor.elevList;
    azimList = indoor.azimList;
    gainTot = indoor.gainTot;
else
    error('ERROR: Wrong environment. use outdoor or indoor')
end


[X,Y] = meshgrid(elevList,azimList);  % Orientation space (exhaustive)

figure;  hold on
fig_idx = get(gcf,'Number');

idx = 1;
hp = [];
hg = [];
for expID = expIDList
    for antID = antIDList
        % get color ID
        colId = mod(((expID-1)*length(antIDList)) + antID,length(colOrder)) + 1;
        % Channel gain (exhaustive)
        Z = gainTot(:,:,antID,expID);
        % calculate the sample variogram
        a = reshape(X,[size(X,1)*size(X,2),1]);
        b = reshape(Y,[size(Y,1)*size(Y,2),1]);
        c = reshape(Z,[size(Z,1)*size(Z,2),1]);
        v = variogram([a b],c,'maxdist',150,'plotit',false);
        [~,~,~,vstruct] = variogramfit(v.distance,v.val,[],[],[],'model','gaussian','plotit',false);
        figure(fig_idx);
        p = plot(vstruct.h,vstruct.gammahat,'color',colOrder(colId,:),'lineWidth',2);
        hp1 = plot(v.distance,v.val,'.','lineStyle','none','lineWidth',2,'MarkerSize',15,'color',colOrder(colId,:));
        hp = [hp hp1];
        hg = [hg p];
        myString{idx} = sprintf('TX %d - RX %d',expID,antID);
        myString{idx+length(antIDList)*length(expIDList)} = sprintf('TX %d - RX %d (Gaussian Fit)',expID,antID);
        idx = idx + 1;
    end
end
legend([hp hg],myString,'location','East Outside');
ylabel('\gamma(h)','FontSize',12);
xlabel('Angular distance','FontSize',12);
title(sprintf('%s - Gaussian fit using BLUE',environment),'FontSize',12);
pos = get(gcf, 'Position');
set(gcf, 'Position', [pos(1) pos(2) 570 306]);
grid minor;