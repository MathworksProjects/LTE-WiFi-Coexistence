%========================= IMPORTANT NOTE ================================%
% ** Author: Carlos Bocanegra Guerra, PhD student
% ** email: bocanegrac@coe.neu.edu
% ** Northeastern University, Boston MA
%
% This code is the runnable version of the code "LTE_WiFi_TxChRx.m", which
% evaluates the throughput, BER and Frame Detection Errors for
% LTE and WiFi standard when coexisting using Almost Blank Subframes (ABS).
% Here, the Tx parameters on the WiFi and LTE side as well as the channel
% conditions can be configured. Giving a certain location in the XY plane
% of the BS and the AP, the WiFi devices are located in a circular sector
% around the AP and the performance is evaluated for every location. For a
% rectangular location of the devices, please run the code from
% "LTE_WiFi_TxChRx_Runnable_rectangular.m".
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

clear; close all;

ABSConfig = 1;              % ABS Configuration

mcs_range = [0, 1, 2, 3];
BER0 = []; BER_glob = [];
not_detected0 = []; not_detected_glob = [];
Prx_WiFi0 = []; Prx_WiFi_glob = [];
Prx_LTE0 = []; Prx_LTE_glob = [];
SINR_WiFi0 = []; SINR_WiFi_glob = [];
PER_WiFi0 = []; PER_WiFi_glob = [];
D_WiFi_glob = [];
D_LTE_glob = [];
BER0_new = []; BER_new_glob = [];

posNodes = [1:79; zeros(1,79)]';

D_WiFi = posNodes(:, 1);
D_LTE  = 80-posNodes(:, 1);

for MCS = mcs_range
        [PSDULength] = wifi_config(MCS);             % Length of the PSDU Field
        tic
        parfor ind = 1:size(D_WiFi,1)
            distance_WiFi = D_WiFi(ind);
            distance_LTE  = D_LTE(ind);
            rng(ind);
            [ber, notDetected, per_Wifi, prx_Wifi, prx_LTE, sinr_Wifi, ber_new] = ...
                LTE_WiFi_TxChRx(distance_WiFi, distance_LTE, ABSConfig, MCS);
            BER0 = [BER0 ber];
            not_detected0 = [not_detected0 notDetected];
            PER_WiFi0 = [PER_WiFi0 per_Wifi];
            Prx_WiFi0 = [Prx_WiFi0 prx_Wifi];
            Prx_LTE0 = [Prx_LTE0 prx_LTE];
            SINR_WiFi0 = [SINR_WiFi0 sinr_Wifi];
            BER0_new = [BER0_new ber_new];
        end
        toc
        D_WiFi_glob = [D_WiFi_glob D_WiFi];                                           %#ok<AGROW>
        D_LTE_glob = [D_LTE_glob D_LTE];                                              %#ok<AGROW>
        BER_glob = [BER_glob BER0]; BER0 = [];                                        %#ok<AGROW>
        not_detected_glob = [not_detected_glob not_detected0]; not_detected0 = [];    %#ok<AGROW>
        Prx_WiFi_glob = [Prx_WiFi_glob Prx_WiFi0]; Prx_WiFi0 = [];                    %#ok<AGROW>
        Prx_LTE_glob = [Prx_LTE_glob Prx_LTE0]; Prx_LTE0 = [];                        %#ok<AGROW>
        SINR_WiFi_glob = [SINR_WiFi_glob SINR_WiFi0]; SINR_WiFi0 = [];                %#ok<AGROW>
        PER_WiFi_glob = [PER_WiFi_glob PER_WiFi0]; PER_WiFi0 = [];                    %#ok<AGROW>
        BER_new_glob = [BER_new_glob BER0_new]; BER0_new = [];                        %#ok<AGROW>
end

save('results_MCS0123_ALL.mat');