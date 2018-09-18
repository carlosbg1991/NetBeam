clear all; close all; clc;

nTxAntennas = 2;

% lastFeedback1 = rand(1,nTxAntennas) + 1j*rand(1,nTxAntennas);
% fid1 = fopen('weights_tx2.bin','wb');
% fwrite(fid1,[real(lastFeedback2) imag(lastFeedback2)],'double');
% fclose(fid1);

maxIter = 1e10;
for iter = 1:maxIter
    tic;
    fid2 = fopen('weights_tx2.bin','rb');
    fseek(fid2,0,'bof');
    channelEst1 = fread(fid2,nTxAntennas*2,'double');
    toc
    fprintf('%.1f\n',channelEst1)
    pause(0.1);
    fclose(fid2);
end