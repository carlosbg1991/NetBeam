% Configure transmitter
radioConfig.platform1            = 'X310';
radioConfig.IPAddr1              = '192.168.10.2';
radioConfig.MasterClockRate1     = 200e6;
radioConfig.Interpolationfactor1 = 500;
radioConfig.ChannelMapping       = 1;
radioConfig.CenterFrequency      = 900e6;
radioConfig.Gain                 = 'External';
radioConfig.ClockSource          = 'External';
radioConfig.PPSSource            = 'External';
% Create transmitter SO
transmitter = comm.SDRuTransmitter('Platform',radioConfig.Platform,...
                                   'IPAddress',radioConfig.IPAddrs,...
                                   'MasterClockRate',radioConfig.MasterClockRate1,...
                                   'InterpolationFactor',radioConfig.Interpolationfactor1,...
                                   'ChannelMapping',radioConfig.ChannelMapping,...
                                   'CenterFrequency',radioConfig.CenterFrequency,...
                                   'Gain',radioConfig.Gain,...
                                   'ClockSource',radioConfig.ClockSource,...
                                   'PPSSource',radioConfig.PPSSource);
% Construct training signals - Golay Sequence 1 (do not modify)
goldSeq = comm.GoldSequence;
goldSeq.FirstPolynomial = [11 2 0];
goldSeq.SecondPolynomial = [11 8 5 2 0];
goldSeq.SamplesPerFrame = 2047;
goldSeq.Index = 15;
goldSeq.FirstInitialConditions = [0 1 0 1 1 0 1 0 0 1 1];
goldSeq.SecondInitialConditions = [0 0 0 1 1 0 1 0 1 1 0];
goldSequence1 = goldSeq();
txfilter = comm.RaisedCosineTransmitFilter;
trainingSig = txfilter([goldSequence1 ; zeros(10,1)]); % Pad extra zeros at the end

%Main loop
for i = 1:30000
    txSig = [trainingSig; zeros(400,4)] * 0.2;
    transmitter(txSig);
end

release(transmitter);

%% Appendix
% This example uses the following helper functions:
%
% * <matlab:edit('helperMUBeamformAllocateRadios.m') helperMUBeamformAllocateRadios.m>
% * <matlab:edit('helperMUBeamformInitGoldSeq.m') helperMUBeamformInitGoldSeq.m>
% * <matlab:edit('helperMUBeamformRx1.m') helperMUBeamformRx1.m>
% * <matlab:edit('helperMUBeamformRx2.m') helperMUBeamformRx2.m>
% * <matlab:edit('helperMUBeamformEstimateChannel.m') helperMUBeamformEstimateChannel.m>

displayEndOfDemoMessage(mfilename) 

%% Copyright Notice
% Universal Software Radio Peripheral(R) and USRP(R) are trademarks of
% National Instruments Corp.