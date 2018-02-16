% 802.11 Channel and Noise
function [ rxWaveform_awgn, Prx_WiFi ] = wiFi_channel(txWaveform, distance, Pnoise, fsWiFi, config, conf, delayProfile, seed)
    %% Fading Channel
    if(strcmp(config,'NonHT'))
        % Channel Model;
        switch(delayProfile)        % Delay Spread selection
            case 'Model-A'; rms = 50e-9;
            case 'Model-B'; rms = 100e-9;
            case 'Model-C'; rms = 150e-9;
            case 'Model-D'; rms = 140e-9;
            case 'Model-E'; rms = 250e-9;
            case 'Model-F'; rms = 150e-9;
        end
        fd = 3;
        channel = stdchan(1/fsWiFi,fd,'802.11a',rms);
        rxWaveform = filter(channel,txWaveform);
    elseif(strcmp(config,'HT'))
        % Create and configure the channel
        tgnChannel = wlanTGnChannel;
        tgnChannel.DelayProfile = delayProfile;
        tgnChannel.NumTransmitAntennas = 1;
        tgnChannel.NumReceiveAntennas = 1;
        tgnChannel.TransmitReceiveDistance = distance; % Distance in meters
        tgnChannel.LargeScaleFadingEffect = 'None';
        tgnChannel.SampleRate = fsWiFi;
        rxWaveform = step(tgnChannel,txWaveform);
        reset(tgnChannel);
    elseif(strcmp(config,'VHT'))
        % Create and configure the channel
        tgacChannel = wlanTGacChannel;
        tgacChannel.DelayProfile = delayProfile;
        tgacChannel.NumTransmitAntennas = 1;
        tgacChannel.NumReceiveAntennas = 1;
        tgacChannel.TransmitReceiveDistance = distance; % Distance in meters
        tgacChannel.LargeScaleFadingEffect = 'None';
        tgacChannel.SampleRate = fsWiFi;
        rxWaveform = step(tgacChannel,txWaveform);
        reset(tgacChannel);
    end

    %% Path Loss - ITU Model for indoor attenuation
    Pf = 0;     % Floor Penetration loss factor for a regular office in the 5GHz band
                % When transmitting and receiving in the same floor.
    N = 31;     % Distance Power Loss Coeficient for an office Area in the 5GHz band
    pathloss = -(20*log10(5.2e3) - 28 + N*log10(distance) + Pf);
    pathloss = 10^(pathloss/10);

    [lWiFi] = wifi_slength(config, conf, 0, helperSampleRate(conf));
    P_WiFi  = coex_power(rxWaveform(1:lWiFi,1));
    Prx_WiFi = pathloss*P_WiFi;
    rxWaveform = sqrt(Prx_WiFi/P_WiFi).*rxWaveform;
    
    %% Additive Noise
    awgnChannel = comm.AWGNChannel;
    awgnChannel.NoiseMethod = 'Variance';
    awgnChannel.VarianceSource = 'Input port';
    nVar = 10^((Pnoise-30)/10);
    rxWaveform_awgn = step(awgnChannel, rxWaveform, nVar);
end