% 802.11 Transmitter
function [txWaveform,txPSDU,conf] = wiFi_tx(config, chBandwidth, modulation, nAntennas, spatialStreams,...
                                            MCS, PSDULength, nWiFiFrames, wiFiFrameStride, Ptx_WiFi)
    if(strcmp(config,'NonHT'));
        conf = wlanNonHTConfig;
        conf.Modulation = modulation;
        conf.PSDULength = PSDULength;
    elseif(strcmp(config,'HT')); 
        conf = wlanHTConfig;
        conf.PSDULength = PSDULength;
        conf.NumSpaceTimeStreams = spatialStreams;
        conf.SpatialMapping = 'Direct';
        conf.GuardInterval = 'Long';
        conf.RecommendSmoothing = true;
    elseif(strcmp(config,'VHT')); 
        conf = wlanVHTConfig;
        conf.NumSpaceTimeStreams = spatialStreams;
        conf.SpatialMapping = 'Direct';
        conf.GuardInterval = 'Long';
        conf.APEPLength = 3000;
    end
    conf.ChannelBandwidth = chBandwidth;
    conf.MCS = MCS;
    conf.NumTransmitAntennas = nAntennas;

    txPSDU = randi([0 1],conf.PSDULength*8,1); % Generate PSDU data in bits
    
    txWaveform = wlanWaveformGenerator(txPSDU,conf,'NumPackets',nWiFiFrames,...
                    'IdleTime',wiFiFrameStride);

    %% Transmitted Power
	Ptx_WiFi = 10^((Ptx_WiFi-30)/10);
    [lWiFi,~] = wifi_slength(config, conf, 0, helperSampleRate(conf));
    P_WiFi  = coex_power(txWaveform(1:lWiFi,1));
    txWaveform = sqrt(Ptx_WiFi/P_WiFi).*txWaveform;
end