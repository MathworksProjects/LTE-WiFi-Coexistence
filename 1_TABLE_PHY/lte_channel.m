function [ rxWaveform_awgn, lte_pad, Prx_LTE] = lte_channel(txWaveform, distance, Pnoise, fsLTE, sf, channel, seed)
    %% Fadding Channel (Rayleigh)
    channel.InitTime = sf/1000;
    lte_pad = 25;
    rxWaveform = lteFadingChannel(channel,[txWaveform ; zeros(lte_pad,size(txWaveform,2))]);

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
    awgnChannel = comm.AWGNChannel;
    awgnChannel.NoiseMethod = 'Variance';
    awgnChannel.VarianceSource = 'Input port';
    nVar = 10^((Pnoise-30)/10);
    rxWaveform_awgn = step(awgnChannel, rxWaveform, nVar);
end