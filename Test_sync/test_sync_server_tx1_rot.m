clear all; close all; clc;

nTxAntennas = 2;

fid1 = fopen('weights_tx1.bin','rb');

maxIter = 1e10;

oldAngle = 0;
a = arduino('/dev/ttyACM1', 'Uno', 'Libraries', 'Servo');
s2 = servo(a, 'D8', 'MinPulseDuration', 550*10^-6, 'MaxPulseDuration', 2500*10^-6); % azimuth servo
for iter = 1:maxIter
    tic;
    fseek(fid1,0,'bof');
    myData = fread(fid1,nTxAntennas*2 + 2,'double');
    timeCorrection = myData(nTxAntennas*2 + 1);
    angleCorrection = abs(myData(nTxAntennas*2 + 2) - oldAngle);
    toc
    fprintf('DATA: %.1f\n',myData(1:nTxAntennas*2));
    fprintf('Time Correction: %.1f\n',timeCorrection);
    fprintf('Rotation Angle: %.1f\n',angleCorrection);
    writePosition(s2, angleCorrection / 360);
    pause(3);
end

fclose(fid1);