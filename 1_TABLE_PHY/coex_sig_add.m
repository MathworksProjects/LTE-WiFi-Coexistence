function [lte_wifi_waveform, nWiFiFrames,offset] = coex_sig_add(lWiFi, lPadding, lte_rxWaveform_awgn, wifi_rxWaveform_awgn, nWiFiFrames, fsLTE, fsWiFi)
        lLTE = length(lte_rxWaveform_awgn);
        lWiFi_tot = lWiFi + lPadding;
        lWiFi_tot_reshape = (fsLTE/fsWiFi)*lWiFi_tot;

        nWiFiframes_MAX = floor(lLTE/lWiFi_tot_reshape);

        if(nWiFiFrames>nWiFiframes_MAX);   nWiFiFrames = nWiFiframes_MAX;   end

        wifi_rxWaveform_awgn = wifi_rxWaveform_awgn(1:nWiFiFrames*lWiFi_tot);

        rng('shuffle');
        range_off_reshape = lLTE - (lWiFi_tot_reshape*nWiFiFrames);
        range_off = floor((fsWiFi/fsLTE)*range_off_reshape);
        offset = randi(range_off);
        rxWiFi = [zeros(offset,1);wifi_rxWaveform_awgn;zeros(range_off-offset,1)];


        [overlap_reshape,~] = resample(rxWiFi,fsLTE,fsWiFi);   % The sampling frequency
                                                               % difference need to be 
                                                               % compensated if signals
                                                               % are manipulated in the 
                                                               % discrete domain
        lte_wifi_waveform = overlap_reshape + lte_rxWaveform_awgn;
%         lte_wifi_waveform = overlap_reshape;
%         close all
%         coex_plotme(overlap_reshape)
%         hold on
%         coex_plotme(lte_rxWaveform_awgn)
%         hold off
end