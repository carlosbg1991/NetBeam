
maxIter = 1e10;
nTxAntennas1 = 2;
nTxAntennas2 = 0;
nTxAntennas = nTxAntennas1 + nTxAntennas2;

fid = fopen('helperMUBeamformfeedback1.bin','wb');

for iter=1:maxIter
    tic
   
    fseek(fid,0,'bof');
    
    chTx1 = (nTxAntennas*2*(iter-1) : ...
             nTxAntennas*2*(iter-1) + nTxAntennas1*2 - 1);
	timeCorr1 = 0;
    angles1 = [mod(iter,90) mod(iter,180)];
    angles1 = repmat(angles1,1,nTxAntennas1);
    chTx2 = (nTxAntennas*2*(iter-1) + nTxAntennas1*2 : ...
             nTxAntennas*2*iter - 1);
	timeCorr2 = 0;
    angles2 = [mod(iter,90) mod(iter,180)];
    angles2 = repmat(angles2,1,nTxAntennas1);    
    for k=1:nTxAntennas1*2
        fprintf('%.1f  ',chTx1(k));
    end
    fprintf('\n');
    fprintf('%.1f\n',timeCorr1);
    fprintf('%.1f\n',angles1);
    for k=1:nTxAntennas2*2
        fprintf('%.1f  ',chTx2(k));
    end
    fprintf('%.1f\n',timeCorr2);
    fprintf('%.1f\n',angles2);
    
    if mod(iter,10)==0
        fwrite(fid,[chTx1 timeCorr1 angles1 ...
                    chTx2 timeCorr2 angles2 ],'double');
    end
    
    toc
   
    pause(0.1);
    
end

fclose(fidtx1);