clear; close all; clc;

% Load the MAC Variables
run('../Configuration_PHY');

% Locations for Received Power Vector [-30:5:90] dBm
wifi_ap_location = [1.06, 1.53, 2.21, 3.21, 4.65, 6.74, 9.77, 14.16, 20.53, 29.75,...
                    43.12, 62.54, 90.64];
% Number of points per power radius for stable results
np               = [2 2 2 2 2 2 4 4 4 8 8 16 16];
% Distance between the AP and BS in meters
lte_bs_location = [(1:1:20) (21:5:100) (101:100:500)];

mcs_range = conf.mcs_range;
ABSConfig = conf.ABSConfig;

for k = 1:length(wifi_ap_location)
    lte_bs_location = setdiff(lte_bs_location,wifi_ap_location(k));
end

for MCS = mcs_range
    
    % 1. *_glob stores all the curresponding values
    % 2. D_WiFi: distance to WiFi AP (m); D_LTE: distance to LTE BS (m); Prx_WiFi:
    % receive power of WiFi (dBm); Prx_LTE: receive power of LTE (dBm); SINR_WiFi: SINR
    % of WiFi devices (dB); PER_WiFi: Packet Error Rate of WiFi devices
    % (percentage)

    D_WiFi_glob = [];
    D_LTE_glob = [];
    BER_glob = []; BER0 = [];
    not_detected_glob = []; not_detected0 = [];
    Prx_WiFi_glob = []; Prx_WiFi0 = [];
    Prx_LTE_glob = []; Prx_LTE0 = [];
    SINR_WiFi_glob = []; SINR_WiFi0 = [];
    PER_WiFi_glob = []; PER_WiFi0 = [];
    
    [conf.PSDULength] = wifi_config(MCS);  % Length of the PSDU Field
    
    for d_lte = lte_bs_location
        
        fprintf('MCS = %d, LTE BS position (x-axis) = %.3fm\n', MCS, d_lte);
        
        % Place the corresponding devices
        D_LTE = []; D_WiFi = [];
        for pow = 1:length(wifi_ap_location)
            angle = pi/(2*(np(pow)-1));
            for p = 1:np(pow)
                pos_X = wifi_ap_location(pow)*cos((p-1)*angle);
                pos_Y = wifi_ap_location(pow)*sin((p-1)*angle);
                D_LTE = [D_LTE sqrt( (d_lte-pos_X).^2 + (0-pos_Y).^2 )];                   %#ok<*AGROW>
                D_WiFi = [D_WiFi sqrt( (pos_X).^2 + (0-pos_Y).^2 )];                      
            end
        end 
        
        % Start simulation and obtain results
        tic
        parfor ind = 1:size(D_WiFi,2)
            distance_WiFi = D_WiFi(ind);
            distance_LTE  = D_LTE(ind);
            rng(ind);
            [ber, notDetected, per_Wifi, prx_Wifi, prx_LTE, sinr_Wifi] = ...
                LTE_WiFi_TxChRx(conf,distance_WiFi, distance_LTE, ABSConfig, MCS);
            BER0 = [BER0 ber];
            not_detected0 = [not_detected0 notDetected];
            PER_WiFi0 = [PER_WiFi0 per_Wifi];
            Prx_WiFi0 = [Prx_WiFi0 prx_Wifi];
            Prx_LTE0 = [Prx_LTE0 prx_LTE];
            SINR_WiFi0 = [SINR_WiFi0 sinr_Wifi];
            
            fprintf('Distance to WiFi AP = %.3fm, Distance to LTE BS = %.3fm, Prx_WiFi = %.3fdBm, Prx_LTE = %.3fdBm, SINR = %.3fdB, BER = %.3f%%, Not detected = %.3f%%,  PER = %.3f%%\n',...
                D_WiFi(ind), D_LTE(ind), prx_Wifi, prx_LTE, sinr_Wifi, ber, notDetected, per_Wifi);
        end
        toc
        
        D_WiFi_glob = [D_WiFi_glob D_WiFi];                                          
        D_LTE_glob = [D_LTE_glob D_LTE];                                             
        BER_glob = [BER_glob BER0]; BER0 = [];                                       
        not_detected_glob = [not_detected_glob not_detected0]; not_detected0 = [];    
        Prx_WiFi_glob = [Prx_WiFi_glob Prx_WiFi0]; Prx_WiFi0 = [];                   
        Prx_LTE_glob = [Prx_LTE_glob Prx_LTE0]; Prx_LTE0 = [];                       
        SINR_WiFi_glob = [SINR_WiFi_glob SINR_WiFi0]; SINR_WiFi0 = [];                
        PER_WiFi_glob = [PER_WiFi_glob PER_WiFi0]; PER_WiFi0 = [];                  
        
        disp('------------------------------------------------------------------------------------------------------------------------------------');
    end

    if MCS==0
        save('MCS0-ABS1-ALL.mat');
        save('MCS0-ABS1.mat','D_WiFi_glob','D_LTE_glob','Prx_WiFi_glob','Prx_LTE_glob','SINR_WiFi_glob','PER_WiFi_glob');
    end
    if MCS==1
        save('MCS1-ABS1-ALL.mat');
        save('MCS1-ABS1.mat','D_WiFi_glob','D_LTE_glob','Prx_WiFi_glob','Prx_LTE_glob','SINR_WiFi_glob','PER_WiFi_glob');
    end
    if MCS==2
        save('MCS2-ABS1-ALL.mat');
        save('MCS2-ABS1.mat','D_WiFi_glob','D_LTE_glob','Prx_WiFi_glob','Prx_LTE_glob','SINR_WiFi_glob','PER_WiFi_glob');
    end
    if MCS==3
        save('MCS3-ABS1-ALL.mat');
        save('MCS3-ABS1.mat','D_WiFi_glob','D_LTE_glob','Prx_WiFi_glob','Prx_LTE_glob','SINR_WiFi_glob','PER_WiFi_glob');
    end
    
end