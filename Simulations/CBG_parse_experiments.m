%% Clean workspace if called directly by user
if length(dbstack) == 1
    clear all; clear classes; close all; clc;  %#ok
end

%% Total param names
paramList = {'azimList','elevList','chTot','nTxAntennas','appliedElev','appliedAzym'};

%% INDOOR DATA
indoor.fileNameList = {'results-25-Nov-2018_20:15:13.mat',...  % Indoors
                      'results-27-Nov-2018_01:32:43.mat',...
                      'results-29-Nov-2018_04:41:46.mat',...
                      'results-29-Nov-2018_22:45:13.mat',...
                      'results-01-Dec-2018_15:24:20.mat',...
                      ...                                     % PAPER
                      'results-09-Dec-2018_03:50:25.mat',...
                      'results-10-Dec-2018_03:38:01.mat',...
                      'results-11-Dec-2018_02:28:34.mat',...
                      'results-13-Dec-2018_19:29:49.mat',...
                      'results-13-Dec-2018_14:41:42.mat',...
                      'results-13-Dec-2018_02:33:53.mat',...
                      'results-11-Dec-2018_21:41:14.mat',...
                      'results-01-Dec-2018_15:24:20.mat',...
                      'results-12-Dec-2018_15:56:54.mat',...
                      };

indoor.azimLoS = [80,55,135,100;...    % Indoors
                  55,40,110,70;...
                  110,75,135,60;...
                  75,75,75,75;...
                  80,65,115,90;...
                  ...
                  90,90,90,90;...          % Fake - to be replaced
                  90,90,90,90;...          % Fake - to be replaced
                  90,90,90,90;...          % Fake - to be replaced
                  90,90,90,90;...          % Fake - to be replaced
                  90,90,90,90;...          % Fake - to be replaced
                  90,90,90,90;...          % Fake - to be replaced
                  90,90,90,90;...          % Fake - to be replaced
                  90,90,90,90;...          % Fake - to be replaced
                  90,90,90,90;...          % Fake - to be replaced
                 ];

indoor.elevLoS = [90,90,90,90,;...     % Indoors
                 68,68,68,68;...
                 60,60,60,60;...
                 65,65,65,65;...
                 70,67,70,67;...
                 ...
                 45,45,45,45;...           % Fake - to be replaced
                 45,45,45,45;...           % Fake - to be replaced
                 45,45,45,45;...           % Fake - to be replaced
                 45,45,45,45;...           % Fake - to be replaced
                 45,45,45,45;...           % Fake - to be replaced
                 45,45,45,45;...           % Fake - to be replaced
                 45,45,45,45;...           % Fake - to be replaced
                 45,45,45,45;...           % Fake - to be replaced
                 45,45,45,45;...           % Fake - to be replaced
                 ];

%% OUTDOOR DATA
outdoor.fileNameList = {'results-03-Dec-2018_13:12:23.mat',...  % Outdoors
                       'results-03-Dec-2018_14:22:34.mat',...
                       'results-03-Dec-2018_15:32:52.mat',...
                       'results-03-Dec-2018_16:41:21.mat',...
                       'results-03-Dec-2018_17:49:02.mat',...
                       ...                                     % PAPER
                       'results-09-Dec-2018_14:41:54.mat',...
                       'results-09-Dec-2018_16:12:08.mat',...
                       'results-09-Dec-2018_17:18:30.mat',...
                       'results-10-Dec-2018_16:25:43.mat',...
                       'results-09-Dec-2018_19:58:46.mat',...
                       'results-09-Dec-2018_18:56:06.mat',...
                       'results-10-Dec-2018_17:58:02.mat',...
                       'results-03-Dec-2018_16:41:21.mat',...
                       'results-10-Dec-2018_18:58:02.mat',...
                       };

outdoor.azimLoS = [77,55,125,103;...       % Outdoors
                  55,40,104,76;...
                  80,80,80,80;...
                  80,40,140,70;...
                  75,60,120,105;...
                  ...
                  90,90,90,90;...          % Fake - to be replaced
                  90,90,90,90;...          % Fake - to be replaced
                  90,90,90,90;...          % Fake - to be replaced
                  90,90,90,90;...          % Fake - to be replaced
                  90,90,90,90;...          % Fake - to be replaced
                  90,90,90,90;...          % Fake - to be replaced
                  90,90,90,90;...          % Fake - to be replaced
                  90,90,90,90;...          % Fake - to be replaced
                  90,90,90,90;...          % Fake - to be replaced
                  ];

outdoor.elevLoS = [90,90,90,90;...         % Outdoors
                  70,70,70,70;...
                  70,70,70,70;...
                  65,70,70,65;...
                  65,65,65,65;...    
                  ...
                  45,45,45,45;...           % Fake - to be replaced
                  45,45,45,45;...           % Fake - to be replaced
                  45,45,45,45;...           % Fake - to be replaced
                  45,45,45,45;...           % Fake - to be replaced
                  45,45,45,45;...           % Fake - to be replaced
                  45,45,45,45;...           % Fake - to be replaced
                  45,45,45,45;...           % Fake - to be replaced
                  45,45,45,45;...           % Fake - to be replaced
                  45,45,45,45;...           % Fake - to be replaced
                  ];

%% PARSER
indoor = parseData(indoor,paramList,37,19,4);  % INDOORS
outdoor = parseData(outdoor,paramList,19,10,4);  % OUTDOORS

%% SAVER
save('RESULTS','indoor','outdoor','paramList');



                  
function myStruct = parseData(myStruct,paramList,azimLength,elevLength,antennas)

    fileNameList = myStruct.fileNameList;
    myStruct.chTot = zeros(azimLength,elevLength,antennas,length(fileNameList));
    myStruct.gainTot = zeros(azimLength,elevLength,antennas,length(fileNameList));
    
    for fileIdx = 1:length(fileNameList)

        % Retrieve file name
        fileName = fileNameList{fileIdx};
        fprintf('Parsing results from: %s...\n',fileName);
        
        % Store variables in local workspace
        mW = load(fileName,paramList{:});

       % convert into 3D the data space
        chGain = abs(mW.chTot);
        [X,Y] = meshgrid(mW.elevList,mW.azimList);
        dim = size(X,1)*size(X,2);

        % Parse exhaustive space
        Z_prel_gain = zeros(mW.nTxAntennas,length(mW.elevList)*length(mW.azimList));  % intermediate gain variable
        Z_gain = zeros(length(mW.azimList),length(mW.elevList),mW.nTxAntennas);  % final gain variable
        Z_prel_tot = zeros(mW.nTxAntennas,length(mW.elevList)*length(mW.azimList));  % intermediate channel variable
        Z_tot = zeros(length(mW.azimList),length(mW.elevList),mW.nTxAntennas);  % final channel variable
        for id = 1:mW.nTxAntennas
            for t = 1:dim
                elev = X(t);
                azym = Y(t);
                idx_elev = find(mW.appliedElev==elev);
                idx_azym = find(mW.appliedAzym==azym);
                idx = intersect(idx_elev,idx_azym);
                Z_prel_gain(id,t) = mean(chGain(id,idx));
                Z_prel_tot(id,t) = mean(mW.chTot(id,idx));
            end
            Z_gain(:,:,id) = reshape(Z_prel_gain(id,:),[size(X,1),size(X,2)]);
            Z_tot(:,:,id) = reshape(Z_prel_tot(id,:),[size(X,1),size(X,2)]);

            myStruct.gainTot(:,:,id,fileIdx) = Z_gain(:,:,id);  % Total gain
            myStruct.chTot(:,:,id,fileIdx) = Z_tot(:,:,id);  % Total channel
        end
    end
    
    myStruct.elevList = mW.elevList;
    myStruct.azimList = mW.azimList;
end