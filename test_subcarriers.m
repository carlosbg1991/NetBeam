NFFT          = 256;  % Point for the FFT (OFDM modulator)
Ndc           = 6;  % Zero padding for central frequencies
Nedge         = 10;  % Zero padding for the edges (multipath effect)    

M = 64;  k = log2(M);  % 64QAM
bits{1} = randi([0 1],64*k,1);
M = 32;  k = log2(M);  % 32QAM
bits{2} = randi([0 1],64*k,1);
M = 16;  k = log2(M);  % 16QAM
bits{3} = randi([0 1],64*k,1);
M = 8;  k = log2(M);  % 8QAM
bits{4} = randi([0 1],64*k,1);
M = 4;  k = log2(M);  % QPSK
bits{5} = randi([0 1],64*k,1);
M = 2;  k = log2(M);  % BPSK
bits{6} = randi([0 1],64*k,1);

symbols(1,:) = qammod(bits{1},64,'InputType','bit','UnitAveragePower',true).';  % 64QAM
symbols(2,:) = qammod(bits{2},32,'InputType','bit','UnitAveragePower',true).';  % 32QAM
symbols(3,:) = qammod(bits{3},16,'InputType','bit','UnitAveragePower',true).';  % 16QAM
symbols(4,:)= qammod(bits{4},8,'InputType','bit','UnitAveragePower',true).';  % 8QAM
symbols(5,:)= qammod(bits{5},4,'InputType','bit','UnitAveragePower',true).';  % QPSK
symbols(6,:)= qammod(bits{6},2,'InputType','bit','UnitAveragePower',true).';  % BPSK

%% Construct payload 1:
% 64 symbols with IFFT length of 256
% Each symbol uses 8 subcarriers
% Subcarrier 5 uses 64-QAM. Each point of 64-QAM is used once.
% Other subcarriers use QPSK
BlockSize = (NFFT/2) - (Ndc/2) - Nedge;
rep = floor(BlockSize/6);
modOut1 = zeros(NFFT,64);
idxFFT(:,1) = [(Ndc/2)+0*rep+1:(Ndc/2)+1*rep (NFFT/2)+0*rep+1:(NFFT/2)+1*rep];  % 64QAM
idxFFT(:,2) = [(Ndc/2)+1*rep+1:(Ndc/2)+2*rep (NFFT/2)+1*rep+1:(NFFT/2)+2*rep];  % 32QAM
idxFFT(:,3) = [(Ndc/2)+2*rep+1:(Ndc/2)+3*rep (NFFT/2)+2*rep+1:(NFFT/2)+3*rep];  % 16QAM
idxFFT(:,4) = [(Ndc/2)+3*rep+1:(Ndc/2)+4*rep (NFFT/2)+3*rep+1:(NFFT/2)+4*rep];  % 8QAM
idxFFT(:,5) = [(Ndc/2)+4*rep+1:(Ndc/2)+5*rep (NFFT/2)+4*rep+1:(NFFT/2)+5*rep];  % QPSK
idxFFT(:,6) = [(Ndc/2)+5*rep+1:(Ndc/2)+6*rep (NFFT/2)+5*rep+1:(NFFT/2)+6*rep];  % BPSK
modOut1(idxFFT(:,1),:) = repmat(symbols(1,:),2*rep,1);
modOut1(idxFFT(:,2),:) = repmat(symbols(2,:),2*rep,1);
modOut1(idxFFT(:,3),:) = repmat(symbols(3,:),2*rep,1);
modOut1(idxFFT(:,4),:)  = repmat(symbols(4,:),2*rep,1);
modOut1(idxFFT(:,5),:)  = repmat(symbols(5,:),2*rep,1);
modOut1(idxFFT(:,6),:)  = repmat(symbols(6,:),2*rep,1);
% OFDM Modulation and append OFDM blocks to create time frame
payload = reshape(ifft(modOut1),[],1);
% Scale time-domain signal appropriately
payload = payload/max(real(payload))*0.5;

% Set up spectrum analyzer and constellation diagram
spAnalyzer = dsp.SpectrumAnalyzer;
spAnalyzer.SampleRate = 400e3;
spAnalyzer.SpectralAverages = 64;
spAnalyzer.Title = 'Receiver 1';
spAnalyzer.YLimits = [-100 -20];