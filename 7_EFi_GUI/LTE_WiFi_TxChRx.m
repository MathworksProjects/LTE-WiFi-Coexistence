% function [overlap_reshape, lte_rxWaveform_awgn] = LTE_WiFi_TxChRx(distance_LTE, distance_WiFi)
%========================= IMPORTANT NOTE ================================%
% ** Author: Carlos Bocanegra Guerra, PhD student
% ** email: bocanegrac@coe.neu.edu
% ** Northeastern University, Boston MA
%
% When trying to evaluate the LTE and WiFi performance for different BS,
% AP, MS and UE locations, please run this code through
% "LTE_WiFi_TxChRx_runnable.m". Use this script when trying to analyze the 
% WiFi and LTE performance for a certain BS, AP and MS location. To do so, 
% please:
% - Uncomment the parameters at the beginning of the script
% - Uncomment theREPORT section at the end of the code to obtain the 
%   performance metrics.
%
% This code evaluates the throughput, BER and Frame Detection Errors for
% LTE and WiFi standard when coexisting using Almost Blank Subframes (ABS).
% WiFi transmissions are always scheduled in ABS, in which LTE only
% transmits Control Signal. LTE Data is scheduled in non-ABS.
%
% The transmit Powers and the received Noise floor can be configured. For
% Unlicensed Bands, the maximum transmit power is 17-23 dBm and the Noise
% Power varies from -105 dBm to -95 dBm.
%
% The Channel is a Slow Fadding Indoor Channel and different scenarios can
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
% - Access to the LTE and WLAN System toolbox.
%
% Code References:
% [1] DownlinkLTEWaveformModelingExample.m
% [2] DownlinkChannelEstimationEqualizationExample.m
% [3] Channel QualityIndicatorReportingExample.m
% [4] TGn Channel Models, IEEE P802.11 Wireless LANs, May 2004
%=========================================================================%
% 
% clear; close all;
% 
% %========================== CONFIG LTE PARAMETERS ========================%
% nDLRB = 100;              % Number Resource Blocks: 6,15,25,50,75,100
% cqiPeriod = 8;            % CQI Report Periodicity (2,5,10, 16, 20,
%                           % 32, 40, 64, 80, 128 or 160 subfrms)
%                           % TS36.101 Table 9.3.2.1.1-1 establishes 8ms
% cqiOffset = 9;            % Offset deriving from cqi-pmi-ConfigurationIndex
%                           % Make sure that the offset is not included in the
%                           % ABSList
% cqiCorrection = 4;        % CQI selection seem to be higher than need it, 
%                           % introducing errors in the frames. Thus, 
%                           % decreasing the throughput.
% numFrames = 1;            % Number of LTE Frames Generated
% ABSConfig = 1;            % ABS Configuration
% %=========================================================================%
% %========================== CONFIG WiFi PARAMETERS =======================%
% config = 'HT';            % NonHT (802.11a/b/g), HT (802.11n), VHT (802.11ac)
% chBandwidth = 'CBW20';    % CBW20, CBW40, CBW80, CBW160
% modulation = 'OFDM';      % OFDM or DSSS
% nAntennas = 1;            % Number of Transmitting Antennas
% spatialStreams = 1;       % Number of Spatial Streams
% MCS = 0;                  % MCS for 802.11 legacy:
%                           %  from 00 to 07 have 1 spatial stream.
%                           % MCS for 802.11n: 
%                           %  from 00 to 07 have 1 spatial stream.
%                           %  from 08 to 15 have 2 spatial streams.
%                           %  from 16 to 23 have 3 spatial streams.
%                           %  from 24 to 31 have 4 spatial streams.
%                           % MCS for 802.11ac:
%                           %  from 00 to 09 have 1 spatial stream.
% % PSDULength = 100;         % Length of the PSDU Field
% nWiFiFrames = 1;          % Number of frames Being Transmitted
% wiFiFrameStride = 20e-6;  % Interval between frames in seconds
% %-------------------------------------------------------------------------%
% [PSDULength] = wifi_config(MCS);  % Length of the PSDU Field
% %=========================================================================%
% %========================== CONFIG LTE-WiFi ==============================%
% Pnoise = -98;             % Noise Power (dBm)
% Ptx_LTE = 17;             % Transmitted Power for LTE in dBm 
% Ptx_WiFi = 17;            % Transmitted Power for WiFi in dBm
% delayProfile = 'Model-B'; % 802.11n  Tgn Models: A, B, C, D or E
%                           % 802.11ac Tgn Models: A, B, C, D or E
% % ------------------------------------------------------------------------%
% pAP = [0 0];              % Location of the AP
% pLTE = [-80 0];           % Location of the BS
% pDev = [-30 0];           % Location of the WiFi Device
% % pGO = [-20 -40];        % Location of a Group owner. The Group owner acts
%                           % as a relay serving de WiFi Device. The Relay 
%                           % assignation is not covered in this code, but it
%                           % can be used to understand the improvement on
%                           % the performance.
% 
% distance_LTE = sqrt((pDev(1)-pLTE(1)).^2 + (pDev(2)-pLTE(2)).^2);
% distance_WiFi = sqrt((pDev(1)-pAP(1)).^2 + (pDev(2)-pAP(2)).^2);
% % distance_WiFi = sqrt((pDev(1)-pGO(1)).^2 + (pDev(2)-pGO(2)).^2);
% 
% % ------------------------------------------------------------------------%
% % LTE Configuration Parameters
% [ABSList, fsLTE, enbConfig, channel, cec] = lte_config(ABSConfig, nDLRB,...
%                                                 cqiOffset, delayProfile);
% % Subframes to run the experiments over
% if length(ABSList)>1; ABSList_loop = 0:(10*numFrames-1);
% else                  ABSList_loop = (ABSList:10:10*numFrames-1);
% end
% % =========================================================================%

CQIReport = [];                       % Reported CQI values non_ABS subframes
CQIReport_ABS = [];                   % Reported CQI values ABS subframes
CQICalc = [];                         % CQI values calculated in each subframe
CQICalc_ABS = [];                     % CQI values calculated in each subframe

totalCRC = [];                        % crc Received in each subframe
totalTBS = [];                        % TBS used to modulate the PDSCH

cqiBuffer = ones(1,cqiPeriod);        % Buffer from which the eNB selects the next CQI
cqiBuffer_ABS = ones(1,cqiPeriod);    % Buffer from which the eNB selects the next CQI
CQI = cqiBuffer(2);                   % Initial CQI value
offsets = 0; offset = 0;              % Initialise frame offset value

txPSDU_tot = [];                      % Transmitted bits per subframe
rxPSDU_tot = [];                      % Received bits per subframe
BER_tot = [];                         % Bit Error Rate per subframe
e_dct_tot = 0;                        % Number of WiFi frames not detected correctly
pkt_ok_tot = 0;                        % Number of Packets Detected and Decoded Successfully
th_wifi_tot = [];
Prx_WiFi_tot = [];
Prx_LTE_tot = [];

for sf = ABSList_loop
    %% LTE Transmitter
    [lte_txWaveform, subframe, enbConfig, TBS, CQI, cqiPtr] = lte_tx(sf, enbConfig, ABSList, cqiPeriod, cqiBuffer, Ptx_LTE);

    %% LTE Channel
    [ lte_rxWaveform_awgn, lte_pad, Prx_LTE] = lte_channel(lte_txWaveform, distance_LTE, Pnoise, fsLTE, sf,channel, sf);

    if( any(enbConfig.NSubframe==ABSList) )        
        %% 802.11 Transmitter
        [wifi_txWaveform,txPSDU,conf] = wiFi_tx(config, chBandwidth, modulation, nAntennas, spatialStreams,...
                 MCS, PSDULength, nWiFiFrames, wiFiFrameStride, Ptx_WiFi);
        fsWiFi = helperSampleRate(conf);

        %% 802.11 Channel
        [wifi_rxWaveform_awgn, Prx_WiFi] = wiFi_channel(wifi_txWaveform, distance_WiFi, Pnoise, fsWiFi, config, conf, delayProfile, sf);

        %% LTE+802.11 Signal Addition
        [lWiFi,lPadding] = wifi_slength(config, conf, wiFiFrameStride, fsWiFi);
        [lte_wifi_waveform, nWiFiFrames, off_wifi,overlap_reshape,lte_rxWaveform_awgn] = coex_sig_add(lWiFi, ... 
            lPadding, lte_rxWaveform_awgn, wifi_rxWaveform_awgn, nWiFiFrames, ...
            fsLTE, fsWiFi);
    else
        lte_wifi_waveform = lte_rxWaveform_awgn;
    end

    %% LTE Receiver
    [cqiBuffer, cqiBuffer_ABS, CQICalc, CQICalc_ABS, crc, offset, offsets] = ...
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
    [BER, e_dct, pkt_ok, th_wifi] = wiFi_report( conf, txPSDU, rxPSDU, off_wifi_tot, off_wifi_dct_tot, nWiFiFrames);
    BER_tot = [BER_tot BER];                                    %#ok<AGROW>
    th_wifi_tot = [th_wifi_tot th_wifi];                        %#ok<AGROW>
    e_dct_tot = e_dct_tot + e_dct;
    pkt_ok_tot = pkt_ok_tot + pkt_ok;
end

%% LTE Report
[medianCQI, medianCQI_ABS, throughput_glob] = lte_report(totalTBS, totalCRC, CQIReport, CQIReport_ABS, numFrames);

% %% REPORT
% fprintf('------------------- REPORT ---------------------\n');
% fprintf('ABS Configuration: %d\n',ABSConfig);
% fprintf('LTE Average CQI in non ABS subframes: %d\n',ceil(mean(CQIReport)));
% fprintf('LTE Average CQI in ABS subframes: %d\n',ceil(mean(CQIReport_ABS)));
% fprintf('LTE Throughput: %.3f Mbps \n',throughput_glob);
% fprintf('WiFi Packets not detected %.3f %% \n',100*e_dct/(nWiFiFrames*numFrames));
% fprintf('WiFi BER: %.3f %% \n',mean(BER_tot));
% fprintf('------------------------------------------------\n');
% fprintf('Prx(noise): %.1fdBm \n',Pnoise);
% fprintf('Ptx(LTE): %.1fdBm Ptx(WiFi): %.1fdBm \n',Ptx_LTE,Ptx_WiFi);
% fprintf('Prx(LTE): %.1fdBm Prx(WiFi): %.1fdBm \n',10*log10(mean(Prx_LTE_tot))+30,10*log10(mean(Prx_WiFi_tot))+30);
% fprintf('d(LTE): %.1fm d(WiFi): %.1fm \n',distance_LTE,distance_WiFi);
% fprintf('------------------------------------------------\n');
% 
% Pw = 10*log10(mean(Prx_WiFi_tot)) + 30;
% Pl = 10*log10(mean(Prx_LTE_tot)) + 30;
% Pnl = 10^((Pnoise-30)/10);
% 
% fprintf('%d %d %.1f %.1f %.3f %.3f %.3f %.2f %.2f %.2f %.1f\n',ABSConfig, MCS,...
%     distance_WiFi, distance_LTE, mean(BER_tot), 100*e_dct_tot/numFrames, ...
%     100*pkt_ok_tot/(numFrames*nWiFiFrames), Pw, Pl, Pnoise, ...
%     10*log10(mean(Prx_WiFi_tot)/(mean(Prx_LTE_tot) + Pnl)));
% 
% th_w = mean(th_wifi_tot)/1e6;
% fprintf('%d %.1f %.1f %.3f %.3f %.2f %.2f %.2f %.1f \n',MCS, distance_WiFi, ...
%     distance_LTE, Pw, Pl, Pnoise, mean(Prx_WiFi_tot)/(mean(Prx_LTE_tot) + Pnl), ...
%     th_w);
% 
% disp(' ');