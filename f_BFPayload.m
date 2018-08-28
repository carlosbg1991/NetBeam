function [txSig,channelEst_new] = f_BFPayload(numTxAntennas,payload,trainingSig,fileName,channelEst_old)

% Open files for reading only
fid1 = fopen(fileName,'rb');
fseek(fid1,0,'bof'); % Read from the beginning of the file

% The channel between 4 TX antennas and 1 RX antenna is modeled
% by 4 complex gains. This approximation works because the
% signal has very narrow bandwidth (400k samples per second).
channelEst = fread(fid1,numTxAntennas*2,'double');

if length(channelEst) == numTxAntennas*2
    % File content has expected length
    channelEst_new = (   channelEst(1:numTxAntennas) + ...
                   1j*channelEst(numTxAntennas+1:numTxAntennas*2)).';
else
    % Use last feedback
    channelEst_new = channelEst_old;
end

% Rotate payload 1 so that receiver 1 does not need to
% correct the phase of payload 1
beamWeight1 = ones(1,numTxAntennas);  % Does not matter
phaseCorrection1 = beamWeight1 * channelEst_new.';  % Beamforming info
phaseCorrection1 = phaseCorrection1/abs(phaseCorrection1);
beamWeight1 = beamWeight1 ./ phaseCorrection1;

% Beamforming for payload
mod_payload = payload*beamWeight1;

% Send signals to the radios
txSig = [trainingSig; zeros(400,numTxAntennas); mod_payload; zeros(100,numTxAntennas)] * 0.2;



% EOF