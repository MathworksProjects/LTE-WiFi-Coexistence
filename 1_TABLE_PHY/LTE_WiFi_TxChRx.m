function [ber, notDetected, per_Wifi, prx_Wifi, prx_LTE, sinr_Wifi, ber_new] = ...
    LTE_WiFi_TxChRx(conf,distance_WiFi, distance_LTE, ABSConfig, MCS)
%========================= IMPORTANT NOTE ================================%
% ** Author: Carlos Bocanegra Guerra, PhD student
% ** email: bocanegrac@coe.neu.edu
% ** Northeastern University, Boston MA
%
% This code evaluates the throughput, BER and Frame Detection Errors for
% LTE and WiFi standard when coexisting using Almost Blank Subframes (ABS).
% WiFi transmissions are always scheduled in ABS, in which LTE only
% transmits Control Signal. LTE Data is scheduled in non-ABS. Please, call 
% LTE_WiFi_TxChRx_Runnable when evaluating the LTE and Wi-Fi performance 
% for different BS, AP, MS and UE locations.
%
% The Transmit Powers and the received Noise floor are tunnable. For the 
% Unlicensed Band, the maximum transmit power is 17-23 dBm and the Noise
% Power varies from -105 dBm to -95 dBm.
%
% The Channel is a Slow Fading Indoor Channel and different scenarios can
% be configured:
% - Model A for a typical office environment, non-line-of-sight (NLOS)
%   conditions, and 50 ns rms delay spread.
% - Model B for a typical large open space and office environments, NLOS
%   conditions, and 100 ns rms delay spread.
% - Model C for a large open space (indoor and outdoor), NLOS conditions,
%   and 150 ns rms delay spread.
% - Model D, same as model C, line-of-sight (LOS) conditions, and 140 ns
%   rms delay spread (10 dB Ricean K-factor at the first delay).
% - Model E for a typical large open space (indoor and outdoor), NLOS
%   conditions, and 250 ns rms delay spread.
%
% The Path loss model ised is the one proposed by the ITU-R for Indoor
% environments.
%
% The Subframes that are intended to be configured as ABS should be
% configured in lte_config. As a Base configurations, ABSConfig from 0 to 9
% shows 1 ABS subframe which corresponds to the value pick in a subframe.
% For example, ABSConfig = 1 configures subframe 1 as an ABS.
%
% The simulations must be split into two groups:
% (1) The first one is the throughput calculation. We set the PSDULength so
% that for different MCS, the packet length varies and more or less packets
% will be fit within a subframe.
% (2) BER and Frame Detection Error calculation. The PSDULength has to be
% set according to the MCS so that we keep the same Packet Length.
%
% System Requeriments:
% - MATLAB and basic student toolbox set.
% - Access to the LTE and WLAN System toolboxes.
%
% Code References:
% [1] DownlinkLTEWaveformModelingExample.m
% [2] DownlinkChannelEstimationEqualizationExample.m
% [3] Channel QualityIndicatorReportingExample.m
% [4] TGn Channel Models, IEEE P802.11 Wireless LANs, May 2004
%=========================================================================%

% Load Local variables
nDLRB           = conf.nDLRB;
cqiPeriod       = conf.cqiPeriod;
cqiOffset       = conf.cqiOffset;
cqiCorrection   = conf.cqiCorrection;
numFrames       = conf.numFrames;
config          = conf.config;
chBandwidth     = conf.chBandwidth;
modulation      = conf.modulation;
nAntennas       = conf.nAntennas;
spatialStreams  = conf.spatialStreams;
nWiFiFrames     = conf.nWiFiFrames;
wiFiFrameStride = conf.wiFiFrameStride;
PSDULength      = conf.PSDULength;
Pnoise          = conf.Pnoise;
Ptx_LTE         = conf.Ptx_LTE;
Ptx_WiFi        = conf.Ptx_WiFi;
delayProfile    = conf.delayProfile;

% Configure LTE-BS
[ABSList, fsLTE, enbConfig, channel, cec] = lte_config(ABSConfig, nDLRB,...
    cqiOffset, delayProfile);

% Subframes over which to run the experiments
if length(ABSList)>1; ABSList_loop = 0:(10*numFrames-1);
else                  ABSList_loop = (ABSList:10:10*numFrames-1);
end
% =========================================================================%

CQIReport = [];                       % Reported CQI values non_ABS subframes
CQIReport_ABS = [];                   % Reported CQI values ABS subframes
CQICalc = [];                         % CQI values calculated in each subframe
CQICalc_ABS = [];                     % CQI values calculated in each subframe

totalCRC = [];                        % crc Received in each subframe
totalTBS = [];                        % TBS used to modulate the PDSCH

cqiBuffer = ones(1,cqiPeriod);        % Buffer from which the eNB selects the next CQI
cqiBuffer_ABS = ones(1,cqiPeriod);    % Buffer from which the eNB selects the next CQI
offsets = 0;                          % Initialise frame offset value

BER_tot = [];                         % Bit Error Rate per subframe
e_dct_tot = 0;                        % Number of WiFi frames not detected correctly
pkt_ok_tot = 0;                       % Number of Packets Detected and Decoded Successfully
th_wifi_tot = [];                     % Throughput in Kbps
Prx_WiFi_tot = [];                    % Total Wi-Fi Received Power in dBm
Prx_LTE_tot = [];                     % Total LTE Received Power in dBm
BER_new_tot = [];

for sf = ABSList_loop
    %% LTE Transmitter
    [lte_txWaveform, ~, enbConfig, TBS, CQI, cqiPtr] = lte_tx(sf, enbConfig, ABSList, cqiPeriod, cqiBuffer, Ptx_LTE);
    
    %% LTE Channel
    [ lte_rxWaveform_awgn, ~, Prx_LTE] = lte_channel(lte_txWaveform, distance_LTE, Pnoise, fsLTE, sf,channel, sf);
    
    if( any(enbConfig.NSubframe==ABSList) )
        %% 802.11 Transmitter
        [wifi_txWaveform,txPSDU,conf] = wiFi_tx(config, chBandwidth, modulation, nAntennas, spatialStreams,...
            MCS, PSDULength, nWiFiFrames, wiFiFrameStride, Ptx_WiFi);
        fsWiFi = helperSampleRate(conf);
        
        %% 802.11 Channel
        [wifi_rxWaveform_awgn, Prx_WiFi] = wiFi_channel(wifi_txWaveform, distance_WiFi, Pnoise, fsWiFi, config, conf, delayProfile, sf);
        
        %% LTE+802.11 Signal Addition
        [lWiFi,lPadding] = wifi_slength(config, conf, wiFiFrameStride, fsWiFi);
        [lte_wifi_waveform, nWiFiFrames, off_wifi] = coex_sig_add(lWiFi, ...
            lPadding, lte_rxWaveform_awgn, wifi_rxWaveform_awgn, nWiFiFrames, ...
            fsLTE, fsWiFi);
    else
        lte_wifi_waveform = lte_rxWaveform_awgn;
    end
    
    %% LTE Receiver
    [cqiBuffer, cqiBuffer_ABS, CQICalc, CQICalc_ABS, crc, ~, offsets] = ...
        lte_rx(sf, enbConfig, ABSList, lte_rxWaveform_awgn, TBS, cqiBuffer, cqiBuffer_ABS, ...
        CQICalc, CQICalc_ABS, cqiCorrection, cqiPeriod, cqiOffset, cqiPtr, cec, offsets);
    totalCRC = [totalCRC crc];                                  %#ok<AGROW>
    totalTBS = [totalTBS TBS];                                  %#ok<AGROW>
    CQIReport = [CQIReport CQI];                                %#ok<AGROW>
    CQIReport_ABS  = [CQIReport_ABS cqiBuffer_ABS(cqiPtr+1)];   %#ok<AGROW>
    Prx_WiFi_tot = [Prx_WiFi_tot Prx_WiFi];                     %#ok<AGROW>
    Prx_LTE_tot = [Prx_LTE_tot Prx_LTE];                        %#ok<AGROW>
    
    %% 802.11 Receiver
    if( any(enbConfig.NSubframe==ABSList) )
        [rxWiFi,~] = resample(lte_wifi_waveform,fsWiFi,fsLTE);
        rxPSDU = [];
        off_wifi_dct_tot = [];
        off_wifi_tot = [];
        for fr = 1:nWiFiFrames
            [rxPSDU1, off_wifi_dct1] = wiFi_rx(rxWiFi, config, conf, off_wifi);
            rxPSDU = [rxPSDU rxPSDU1];                           %#ok<AGROW>
            off_wifi_dct_tot = [off_wifi_dct_tot off_wifi_dct1]; %#ok<AGROW>
            off_wifi_tot = [off_wifi_tot ; off_wifi];            %#ok<AGROW>
            rxWiFi(1:off_wifi+lWiFi+lPadding,:) = [];
            if(fr == 1); off_wifi = 0; end
        end
    end
    [BER, e_dct, pkt_ok, th_wifi, BER_new] = wiFi_report( conf, txPSDU, rxPSDU, off_wifi_tot, off_wifi_dct_tot, nWiFiFrames);
    BER_tot = [BER_tot BER];                                    %#ok<AGROW>
    th_wifi_tot = [th_wifi_tot th_wifi];                        %#ok<AGROW>
    e_dct_tot = e_dct_tot + e_dct;
    pkt_ok_tot = pkt_ok_tot + pkt_ok;
    BER_new_tot = [BER_new_tot BER_new];                        %#ok<AGROW>
end

%% LTE Report
% [medianCQI, medianCQI_ABS, throughput_glob] = lte_report(totalTBS, totalCRC, CQIReport, CQIReport_ABS, numFrames);

%% RESULTS
ber = mean(BER_tot);
notDetected = 100*e_dct_tot/(numFrames*nWiFiFrames);
per_Wifi = 100*(1-(pkt_ok_tot/(numFrames*nWiFiFrames)));
prx_Wifi = 10*log10(mean(Prx_WiFi_tot)) + 30;
prx_LTE = (10*log10(mean(Prx_LTE_tot)) + 30);
sinr_Wifi = 10*log10(mean(Prx_WiFi_tot)/(mean(Prx_LTE_tot)+(10^((Pnoise-30)/10))));
ber_new = mean(BER_new_tot(~isnan(BER_new_tot)));