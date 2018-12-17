clear all; close all; clc;

nTxAntennas = 2;

fid1 = fopen('weights_tx1.bin','rb');

maxIter = 1e10;
for iter = 1:maxIter
    tic;
    fseek(fid1,0,'bof');
    myData = fread(fid1,nTxAntennas*2 + 1,'double');
    toc
    fprintf('%.1f\n',myData)
    pause(0.1);
end

fclose(fid1);