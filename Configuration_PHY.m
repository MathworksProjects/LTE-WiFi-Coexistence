%========================== CONFIG LTE PARAMETERS ========================%
conf.nDLRB = 100;              % Number Resource Blocks: 6,15,25,50,75,100
conf.cqiPeriod = 8;            % CQI Report Periodicity (2,5,10, 16, 20,
                               % 32, 40, 64, 80, 128 or 160 subfrms)
                               % TS36.101 Table 9.3.2.1.1-1 establishes 8ms
conf.cqiOffset = 9;            % Offset deriving from cqi-pmi-ConfigurationIndex
                               % Make sure that the offset is not included in the
                               % ABSList
conf.cqiCorrection = 4;        % CQI selection seem to be higher than need it,
                               % introducing errors in the frames. Thus,
                               % decreasing the throughput.
conf.numFrames = 200;          % Number of LTE Frames Generated
%=========================================================================%
%========================== CONFIG WiFi PARAMETERS =======================%
conf.config = 'HT';            % NonHT (802.11a/b/g), HT (802.11n), VHT (802.11ac)
conf.chBandwidth = 'CBW20';    % CBW20, CBW40, CBW80, CBW160
conf.modulation = 'OFDM';      % OFDM or DSSS
conf.nAntennas = 1;            % Number of Transmitting Antennas
conf.spatialStreams = 1;       % Number of Spatial Streams
conf.nWiFiFrames = 1;          % Number of frames Being Transmitted
conf.wiFiFrameStride = 20e-6;  % Interval between frames in seconds
%=========================================================================%
%========================== CONFIG LTE-WiFi ==============================%
conf.Pnoise = -98;             % Noise Power (dBm)
conf.Ptx_LTE = 17;             % Transmitted Power for LTE in dBm
conf.Ptx_WiFi = 17;            % Transmitted Power for WiFi in dBm
conf.delayProfile = 'Model-B'; % 802.11n  Tgn Models: A, B, C, D or E
                          % 802.11ac Tgn Models: A, B, C, D or E
% ------------------------------------------------------------------------%

conf.ABSConfig = 1;            % ABS Configuration
conf.mcs_range = [0, 1, 2, 3]; % MCS Range
                          	   % MCS Configuration for 802.11n:
                          	   %  from 00 to 07 have 1 spatial stream.
                               %  from 08 to 15 have 2 spatial streams.
                               %  from 16 to 23 have 3 spatial streams.
                               %  from 24 to 31 have 4 spatial streams.
                               % MCS Configuration for 802.11ac:
                               %  from 00 to 09 have 1 spatial stream.