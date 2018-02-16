function [ sig_length, padding_length ] = wifi_slength( config, conf, wiFiFrameStride, Rs )
pre = [wlanLSTF(conf);wlanLLTF(conf);wlanLSIG(conf)];
    preamble = length(pre);        
    if(strcmp(config,'HT'))
        pre1 = [wlanHTSIG(conf);wlanHTSTF(conf);wlanHTLTF(conf)];
        preamble = preamble + length(pre1);
    elseif(strcmp(config,'VHT'))
        pre1 = [wlanVHTSIGA(conf);wlanVHTSTF(conf);wlanVHTLTF(conf);wlanVHTSIGB(conf)];
        preamble = preamble + length(pre1);
    end

    % Calculate number of OFDM symbols
    if(strcmp(config,'NonHT'))
        [~, ~, ~, ndbps, ~, ~] = getValues(config, conf);
        numTailBits = 6; numServiceBits = 16;
        numDataSymbols = ceil((8*conf.PSDULength + numServiceBits + numTailBits)/ndbps);
    elseif(strcmp(config,'HT'))
        [~, ~, ~, ndbps, nes, nss] = getValues(config, conf);
        numSTS = conf.NumSpaceTimeStreams;
        STBC = numSTS - nss;
        mSTBC = 1 + (STBC~=0);
        numTailBits = 6;
        numDataSymbols = mSTBC * ceil((8*conf.PSDULength + 16 + ...
            numTailBits*nes)/(mSTBC*ndbps));
    else
        APEPLen  = repmat(conf.APEPLength, 1, conf.NumUsers/length(conf.APEPLength));
        [~, ~, ~, ndbps, nes, ~] = getValues(config, conf);
        numTailBits = 6;
        mSTBC = (conf.NumUsers == 1) * (conf.STBC ~= 0) + 1;
        numDataSymbols = max(mSTBC * ceil((8*APEPLen + 16 + ...
        numTailBits * nes)./(mSTBC * ndbps)));
    end

    sig_length = preamble + numDataSymbols*80;
    padding_length = wiFiFrameStride*Rs;
end

function [rate, modulation, nsd, ndbps, nes, nss] = getValues(config, conf)
    if(strcmp(config,'HT'));	MCS = mod(conf.MCS,8);
    else                        MCS = conf.MCS;
    end
    rate_NonHT  = [1/2 3/4 1/2 3/4 1/2 3/4 2/3 3/4  0   0 ];
    rate_HT     = [1/2 1/2 3/4 1/2 3/4 2/3 3/4 5/6  0   0 ];
    rate_VHT    = [1/2 1/2 3/4 1/2 3/4 2/3 3/4 5/6 3/4 5/6];
    modul_NonHT = [1 1 2 2 4 4 6 6 0 0];
    modul_HT    = [1 2 2 4 4 6 6 6 0 0];
    modul_VHT   = [1 1 2 2 4 4 6 6 8 8];

    mapRate  = containers.Map({'NonHT', 'HT', 'VHT'},{rate_NonHT, rate_HT, rate_VHT});
    mapModul = containers.Map({'NonHT', 'HT', 'VHT'},{modul_NonHT, modul_HT, modul_VHT});

    rates = mapRate(config);
    rate = rates(MCS+1);
    modulations = mapModul(config);
    modulation = modulations(MCS+1);
    
    if(strcmp(config,'NonHT'))
        nsd = 48;
        ncbps = nsd * modulation;
        ndbps = ncbps * rate;
        nes = 0; nss = 0;
    elseif(strcmp(config,'HT'))
        switch conf.ChannelBandwidth
            case 'CBW20'
                nsd = 52;
            case 'CBW40'
                nsd = 108;
        end
        nss = floor(conf.MCS/8)+1;
        ncbps = nsd * modulation * nss;
        ndbps = ncbps * rate;
        nes = ceil(ndbps/(4*300));
    else % case VHT
        switch conf.ChannelBandwidth
            case 'CBW20'
                nsd = 52;
            case 'CBW40'
                nsd = 108;
            case {'CBW80', 'CBW80+80'}
                nsd = 234;
            otherwise % case 'CBW160'
                nsd = 468;
        end
        nss = conf.NumSpaceTimeStreams / (((conf.NumUsers == 1) && conf.STBC) + 1);
        ncbps = nsd * modulation * nss;
        ndbps = ncbps * rate;
        
        % Handle exceptions to Nes generic rule - Table 7.13 [2].
        % For each case listed, work off the Ndbps value and create a look-up
        % table for the Nes value.
        % Only 9360 has a valid value from the generic rule also, 
        % all others are exceptions
        NdbpsVec = [2457 8190 9828 9360 14040 9828 16380 19656 21840 14976 22464];
        expNes =   [   3    6    6    6     8    6     9    12    12     8    12];

        exceptIdx = find(ndbps == NdbpsVec);
        if ~isempty(exceptIdx)
            if (ndbps == 9360) && (nss == 5) % One valid case for 160, 80+80
                nes = 5;
            else  % Two exception cases
                nes = expNes(exceptIdx(1));
            end
        else  % Generic rule: 3.6*600 - for a net 600Mbps per encoder
            nes = ceil(ndbps/2160);
        end
    end
end