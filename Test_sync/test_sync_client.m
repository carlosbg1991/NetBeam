
maxIter = 1e10;
nTxAntennas1 = 0;
nTxAntennas2 = 2;
nTxAntennas = nTxAntennas1 + nTxAntennas2;

fid = fopen('helperMUBeamformfeedback1.bin','wb');

for iter=1:maxIter
    tic
   
    fseek(fid,0,'bof');
    
    chTx1 = (nTxAntennas*2*(iter-1) : ...
             nTxAntennas*2*(iter-1) + nTxAntennas1*2 - 1);
    chTx2 = (nTxAntennas*2*(iter-1) + nTxAntennas1*2 : ...
             nTxAntennas*2*iter - 1);
    
    for k=1:nTxAntennas1*2
        fprintf('%.1f  ',chTx1(k));
    end
    fprintf('\n');
    for k=1:nTxAntennas2*2
        fprintf('%.1f  ',chTx2(k));
    end
    fprintf('\n');
    
    fwrite(fid,[chTx1 chTx2],'double');
    
    toc
   
    pause(0.1);
    
end

fclose(fidtx1);