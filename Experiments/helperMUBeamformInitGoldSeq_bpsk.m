function refBPSK = helperMUBeamformInitGoldSeq_bpsk
% Helper function for MultiUserBeamformingExample, helperMUBeamformRx1 and
% helperMUBeamformRx2
%
% Copyright 2016 The MathWorks, Inc.

%% Construct Gold sequence 1
goldSeq = comm.GoldSequence;
goldSeq.FirstPolynomial = [11 2 0];
goldSeq.SecondPolynomial = [11 8 5 2 0];
goldSeq.SamplesPerFrame = 2047;
goldSeq.Index = 15;
goldSeq.FirstInitialConditions = [0 1 0 1 1 0 1 0 0 1 1];
goldSeq.SecondInitialConditions = [0 0 0 1 1 0 1 0 1 1 0];

goldSequence1 = goldSeq();

%% Construct Gold sequence 2
goldSeq = comm.GoldSequence;
goldSeq.FirstPolynomial = [11 2 0];
goldSeq.SecondPolynomial = [11 8 5 2 0];
goldSeq.SamplesPerFrame = 2047;
goldSeq.Index = 26;
% goldSeq.Index = 23;
goldSeq.FirstInitialConditions = [1 1 0 1 0 1 0 0 0 0 0];
goldSeq.SecondInitialConditions = [0 1 0 1 1 1 1 1 1 1 0];

goldSequence2 = goldSeq();

%% Construct Gold sequence 3
goldSeq = comm.GoldSequence;
goldSeq.FirstPolynomial = [11 2 0];
goldSeq.SecondPolynomial = [11 8 5 2 0];
goldSeq.SamplesPerFrame = 2047;
goldSeq.Index = 37;
goldSeq.FirstInitialConditions = [0 0 1 0 0 1 1 0 1 1 1];
goldSeq.SecondInitialConditions = [0 1 1 1 1 0 1 0 0 1 1];

goldSequence3 = goldSeq();

%% Construct Gold sequence 4
goldSeq = comm.GoldSequence;
goldSeq.FirstPolynomial = [11 2 0];
goldSeq.SecondPolynomial = [11 8 5 2 0];
goldSeq.SamplesPerFrame = 2047;
goldSeq.Index = 78;
goldSeq.FirstInitialConditions = [1 0 1 0 1 0 0 0 1 1 0];
goldSeq.SecondInitialConditions = [1 0 0 1 0 0 0 1 0 0 0];

goldSequence4 = goldSeq();

%% Construct Gold sequence 5
goldSeq = comm.GoldSequence;
goldSeq.FirstPolynomial = [11 2 0];
goldSeq.SecondPolynomial = [11 8 5 2 0];
goldSeq.SamplesPerFrame = 2047;
goldSeq.Index = 102;
goldSeq.FirstInitialConditions = [0 0 1 1 0 1 0 1 1 0 1];
goldSeq.SecondInitialConditions = [0 0 0 1 0 1 1 1 0 1 0];

goldSequence5 = goldSeq();

%% Construct Gold sequence 6
goldSeq = comm.GoldSequence;
goldSeq.FirstPolynomial = [11 2 0];
goldSeq.SecondPolynomial = [11 8 5 2 0];
goldSeq.SamplesPerFrame = 2047;
goldSeq.Index = 1;
goldSeq.FirstInitialConditions = [1 1 1 0 0 0 1 0 1 0 1];
goldSeq.SecondInitialConditions = [0 0 1 0 1 1 0 0 1 0 0];

goldSequence6 = goldSeq();

%% Use BPSK
gs1 = (2*goldSequence1-1)*exp(1j*pi/4);
gs2 = (2*goldSequence2-1)*exp(1j*pi/4);
gs3 = (2*goldSequence3-1)*exp(1j*pi/4);
gs4 = (2*goldSequence4-1)*exp(1j*pi/4);
gs5 = (2*goldSequence5-1)*exp(1j*pi/4);
gs6 = (2*goldSequence6-1)*exp(1j*pi/4);

%% Reference signal for cross correlation
% txfilter = comm.RaisedCosineTransmitFilter;
% xRef = txfilter([gs1 gs2 gs3 gs4 gs5 gs6; ...
%                  zeros(10,6)]); % Pad extra zeros

% refBin = [goldSequence1 goldSequence2 goldSequence3 goldSequence4 goldSequence5 goldSequence6];
refBPSK = [gs1 gs2 gs3 gs4];
