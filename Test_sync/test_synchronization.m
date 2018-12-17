clear all; close all; clc;

% lastFeedback1 = rand(1,4) + 1j*rand(1,4);
% fid1 = fopen('mierda.bin','wb');
% fwrite(fid1,[real(lastFeedback1) imag(lastFeedback1)],'double');
% fclose(fid1);

fid1 = fopen('mierda.bin','rb');

maxIter = 1e10;
for iter = 1:maxIter
    tic;
    fseek(fid1,0,'bof');
    channelEst1 = fread(fid1,8,'double');
    toc
    fprintf('%.1f %.1f %.1f %.1f %.1f %.1f %.1f %.1f\n',channelEst1)
    pause(0.1);
end

fclose(fid1);