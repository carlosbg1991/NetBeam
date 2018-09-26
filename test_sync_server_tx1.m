clear all; close all; clc;

nTxAntennas = 1;

% lastFeedback1 = rand(1,nTxAntennas) + 1j*rand(1,nTxAntennas);
% fid1 = fopen('weights_tx1.bin','wb');
% fwrite(fid1,[real(lastFeedback1) imag(lastFeedback1)],'double');
% fclose(fid1);

fid1 = fopen('weights_tx1.bin','rb');

maxIter = 1e10;
for iter = 1:maxIter
    tic;
    fseek(fid1,0,'bof');
    channelEst1 = fread(fid1,nTxAntennas*2,'double');
    toc
    fprintf('%.1f\n',channelEst1)
%     pause(0.1);
    duration = 100;
    java.lang.Thread.sleep(duration*1000);
end

fclose(fid1);