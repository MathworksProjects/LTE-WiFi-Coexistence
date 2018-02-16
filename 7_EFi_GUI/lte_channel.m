function [ rxWaveform_awgn, lte_pad, Prx_LTE] = lte_channel(txWaveform, distance, Pnoise, fsLTE, sf, channel, seed)
    %% Fadding Channel (Rayleigh)
%     lte_pad = 25;    
%     channel = wlanTGnChannel('SampleRate',fsLTE,'LargeScaleFadingEffect',...
%         'None','DelayProfile',delayProfile,'RandomStream','mt19937ar with seed',...
%         'Seed', seed);
%     rxWaveform = step(channel,[txWaveform ; zeros(lte_pad,size(txWaveform,2))]);
    channel.InitTime = sf/1000;
    lte_pad = 25;
    rxWaveform = lteFadingChannel(channel,[txWaveform ; zeros(lte_pad,size(txWaveform,2))]);

    %% Path Loss - TGn Path Loss Model
    % TGn Path Loss Model. Information about this model can be found in
    % IEEE P802.11 Wireless LANs. The model defines a break point distance
    % from which the path loss model switches from Free space to free space
    % with some corrections. The break point distance varies across the
    % Channel Model being used.

%     switch(delayProfile)
%         case 'Model-A'
%             dBP = 5;
%         case 'Model-B'
%             dBP = 5;
%         case 'Model-C'
%             dBP = 5;
%         case 'Model-D'
%             dBP = 10;
%         case 'Model-E'
%             dBP = 20;
%         case 'Model-F'
%             dBP = 30;
%     end
% 
%     if(distance<dBP)
%         pathloss = 147.55 - 20*log10(distance) - 20*log10(5.2e6);
%     else
%         pathloss = 147.55 - 20*log10(dBP) - 20*log10(5.2e6) - 35*log10(distance/dBP);
%     end

    %% Path Loss - ITU Model for indoor attenuation
    Pf = 0;     % Floor Penetration loss factor for a regular office in the 5GHz band
                % When transmitting and receiving in the same floor.
    N = 31;     % Distance Power Loss Coeficient for an office Area in the 5GHz band
    pathloss = -(20*log10(5.2e3) - 28 + N*log10(distance) + Pf);
    pathloss = 10^(pathloss/10);

    sym = 4;    % Reference Signals are sent and can be used to calculate
                % the Received LTE Power. They are located in symbols
                % 0, 4, 7 and 11
    lLTE = length(rxWaveform);
    lte_ini_power = 1+sym*floor((lLTE-lte_pad)/14);
    lte_end_power = (sym+1)*floor((lLTE-lte_pad)/14);

    P_LTE  = coex_power(rxWaveform(lte_ini_power:lte_end_power,1));
    Prx_LTE = pathloss*P_LTE;
	rxWaveform = sqrt(Prx_LTE/P_LTE).*rxWaveform;

	%% Additive Noise
    chAWGN = comm.AWGNChannel('NoiseMethod','Variance','VarianceSource','Input port',...
        'RandomStream','mt19937ar with seed','Seed', seed);
    nVar = 10^((Pnoise-30)/10);
    rxWaveform_awgn = step(chAWGN, rxWaveform, nVar);

%     sa = dsp.SpectrumAnalyzer('SampleRate',fsLTE, ...
%             'ShowLegend',true, ...
%             'Window', 'Rectangular', ...
%             'SpectralAverages',10, ...
%             'ChannelNames',{'Tx','Rx','Rx AWGN'});
%     step(sa,[[txWaveform ; zeros(lte_pad,size(txWaveform,2))] rxWaveform rxWaveform_awgn]);
end