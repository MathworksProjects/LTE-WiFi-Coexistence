% 802.11 Receiver
function [ rxPSDU, pktOff] = wiFi_rx(rxWaveform_awgn, config, conf, pktOffsetOK)
ind = wlanFieldIndices(conf);
PFO = comm.PhaseFrequencyOffset;
PFO.SampleRate = helperSampleRate(conf);
PFO.PhaseOffset = 0;
PFO.FrequencyOffsetSource = 'Input port';

if(strcmp(config,'NonHT'))
    %% DETECT BEGINNING OF WIFI FRAME
    % Packet Detection, Short Pre removal and Coarse Freq Offset
    % Correct detection implies the detection of both, the short and the
    % long preamble. If no preamble is detected, the code returns
    % -999 as an error. The code keeps iterating until it finds the frame
    % or no detection of the short preamble is found.
    rxWaveform_awgn_copy = rxWaveform_awgn;
    off_add = 0;
    found = false;
    while(found == false)
        pktStartIdx = helperPacketDetect(rxWaveform_awgn_copy,conf.ChannelBandwidth);
        if~isempty(pktStartIdx)
            pktOffset = pktStartIdx-1;
            lstf = rxWaveform_awgn_copy(pktOffset+(ind.LSTF(1):ind.LSTF(2)),:);
            coarseFreqOff = wlanCoarseCFOEstimate(lstf,conf.ChannelBandwidth);
            rxWaveform_awgn_corr = step(PFO,rxWaveform_awgn_copy,-coarseFreqOff);
            release(PFO);

            % Long Pre removal and Fine Freq Offset
            nonhtfields = rxWaveform_awgn_corr(pktOffset+(ind.LSTF(1):ind.LSIG(2)),:);
            lltfIdx = helperSymbolTiming(nonhtfields,conf.ChannelBandwidth);
            pktOff = pktOffset+lltfIdx-double(ind.LLTF(1));
            if isempty(lltfIdx) || pktOff<0
                % Short Preamble detection is incorrect. New iteration.
                rxWaveform_awgn_copy(1:300) = [];
                off_add = off_add + 300;
            else
                % Short and Long Preamble are detected
                found = true;
            end
        else
            % Frame not found
            pktOff = -999;
            break;
        end
    end
    
    if(found == true)
        pktOff = pktOff + off_add;
    end
    
    %% DECODE BITS ASSUMING PERFECT SYNCHRONIZATION (KNOWN FRAME START)
    % Coarse and Fine CFO Estimate
    lstf = rxWaveform_awgn(pktOffsetOK+(ind.LSTF(1):ind.LSTF(2)),:);
    coarseFreqOff = wlanCoarseCFOEstimate(lstf,conf.ChannelBandwidth);
    rxWaveform_awgn = step(PFO,rxWaveform_awgn,-coarseFreqOff);
    release(PFO);
    rxWaveform_awgn = rxWaveform_awgn(1+pktOffsetOK:end,:);
    lltf = rxWaveform_awgn(ind.LLTF(1):ind.LLTF(2),:); 
    fineFreqOff = wlanFineCFOEstimate(lltf,conf.ChannelBandwidth);
    rxWaveform_awgn = step(PFO,rxWaveform_awgn,-fineFreqOff);
    release(PFO);

    % Noise Power Estimation and Channel Estimation
    demodLLTF = wlanLLTFDemodulate(lltf,conf);
    nVarNonHT = helperNoiseEstimate(demodLLTF);
    chEst = wlanLLTFChannelEstimate(demodLLTF,conf);

    % Recover Data
    if(ind.NonHTData(2)>length(rxWaveform_awgn))
        rxPSDU = [];
        return;
    end
    rxData = rxWaveform_awgn(ind.NonHTData(1):ind.NonHTData(2),:);
    rxPSDU = wlanNonHTDataRecover(rxData,chEst,nVarNonHT,conf);
elseif(strcmp(config,'HT'))
    %% DETECT BEGINNING OF WIFI FRAME
    % Packet Detection, Short Pre removal and Coarse Freq Offset
    % Correct detection implies the detection of both, the short and the
    % long preamble. If no preamble is detected, the code returns
    % -999 as an error. The code keeps iterating until it finds the frame
    % or no detection of the short preamble is found.
    rxWaveform_awgn_copy = rxWaveform_awgn;
    off_add = 0;
    found = false;
    failCheck = 1;
    while(found == false)
        pktStartIdx = helperPacketDetect(rxWaveform_awgn_copy,conf.ChannelBandwidth);
        if~isempty(pktStartIdx)
            pktOffset = pktStartIdx-1;
            lstf = rxWaveform_awgn_copy(pktOffset+(ind.LSTF(1):ind.LSTF(2)),:);
            coarseFreqOff = wlanCoarseCFOEstimate(lstf,conf.ChannelBandwidth);
%             rxWaveform_awgn_corr = helperFrequencyOffset(rxWaveform_awgn_copy,PFO.SampleRate,-coarseFreqOff);
            rxWaveform_awgn_corr = step(PFO,rxWaveform_awgn_copy,-coarseFreqOff);
            release(PFO);

            % Long Pre removal and Fine Freq Offset
            nonhtfields = rxWaveform_awgn_corr(pktOffset+(ind.LSTF(1):ind.LSIG(2)),:);
            lltfIdx = helperSymbolTiming(nonhtfields,conf.ChannelBandwidth);
            pktOff = pktOffset+lltfIdx-double(ind.LLTF(1));
            if isempty(lltfIdx) || pktOff<0
                % Short Preamble detection is incorrect. New iteration.
                rxWaveform_awgn_copy(1:300) = [];
                off_add = off_add + 300;
            else
                % Coarse and Fine CFO Estimate
                rxWaveform_awgn_corr = rxWaveform_awgn_corr(1+pktOff:end,:);
                lltf = rxWaveform_awgn_corr(ind.LLTF(1):ind.LLTF(2),:); 
                fineFreqOff = wlanFineCFOEstimate(lltf,conf.ChannelBandwidth);
%                 rxWaveform_awgn_corr = helperFrequencyOffset(rxWaveform_awgn_copy,PFO.SampleRate,-fineFreqOff);
                rxWaveform_awgn_corr = step(PFO,rxWaveform_awgn_corr,-fineFreqOff);
                release(PFO);

                % Estimate noise power in HT fields
                lltf = rxWaveform_awgn_corr(ind.LLTF(1):ind.LLTF(2),:);
                demodLLTF = wlanLLTFDemodulate(lltf,conf.ChannelBandwidth);
                nVarHT = helperNoiseEstimate(demodLLTF,conf);

                % Parity check
                chEstNonHT = wlanLLTFChannelEstimate(demodLLTF,conf);
                rxLSIG = rxWaveform_awgn_corr(ind.LSIG(1):ind.LSIG(2),:);
                [~,failCheck,~, ~] = wlanLSIGRecover( ...
                  rxLSIG, chEstNonHT, nVarHT, conf.ChannelBandwidth);

                % Short and Long Preamble are detected
                if(failCheck==0)
                    found = true;
                else
                    % Short Preamble detection is incorrect. New iteration.
                    rxWaveform_awgn_copy(1:300) = [];
                    off_add = off_add + 300;
                end
            end
        else
            % Frame not found
            pktOff = -999;
            break;
        end
    end

    if(found == true)
        pktOff = pktOff + off_add;
    end
    
    %% DECODE BITS ASSUMING PERFECT SYNCHRONIZATION (FRAME START KNOWN)
    % Coarse and Fine CFO Estimate
    lstf = rxWaveform_awgn(pktOffsetOK+(ind.LSTF(1):ind.LSTF(2)),:);
    coarseFreqOff = wlanCoarseCFOEstimate(lstf,conf.ChannelBandwidth);
    rxWaveform_awgn = step(PFO,rxWaveform_awgn,-coarseFreqOff);
    release(PFO);
    rxWaveform_awgn = rxWaveform_awgn(1+pktOffsetOK:end,:);
    lltf = rxWaveform_awgn(ind.LLTF(1):ind.LLTF(2),:); 
    fineFreqOff = wlanFineCFOEstimate(lltf,conf.ChannelBandwidth);
    rxWaveform_awgn = step(PFO,rxWaveform_awgn,-fineFreqOff);
    release(PFO);

    % Estimate noise power in HT fields
    lltf = rxWaveform_awgn(ind.LLTF(1):ind.LLTF(2),:);
    demodLLTF = wlanLLTFDemodulate(lltf,conf.ChannelBandwidth);
    nVarHT = helperNoiseEstimate(demodLLTF,conf);

    % Extract HT-LTF samples from the waveform, demodulate and perform
    % channel estimation
    htltf = rxWaveform_awgn(ind.HTLTF(1):ind.HTLTF(2),:);
    htltfDemod = wlanHTLTFDemodulate(htltf,conf);
    chanEst = wlanHTLTFChannelEstimate(htltfDemod,conf);

    % Recover the transmitted PSDU in HT Data
    % Extract HT Data samples from the waveform and recover the PSDU
    rxData = rxWaveform_awgn(ind.HTData(1):ind.HTData(2),:);
    rxPSDU = wlanHTDataRecover(rxData,chanEst,nVarHT,conf);
elseif(strcmp(config,'VHT'))
    % Packet detect
    pktStartIdx = helperPacketDetect(rxWaveform_awgn,conf.ChannelBandwidth);
    if isempty(pktStartIdx) % If empty no L-STF detected; packet error
        rxPSDU = [];
        return;
    end
    pktOffset = pktStartIdx-1; % Packet offset from start of waveform

    % Extract L-STF and perform coarse frequency offset correction
    lstf = rxWaveform_awgn(pktOffset+(ind.LSTF(1):ind.LSTF(2)),:); 
    coarseFreqOff = wlanCoarseCFOEstimate(lstf,conf.ChannelBandwidth);
    rxWaveform_awgn = step(PFO,rxWaveform_awgn,-coarseFreqOff);
    release(PFO); % Release object for subsequent processing

    % Extract the Non-HT fields and determine start of L-LTF
    nonhtfields = rxWaveform_awgn(pktOffset+(ind.LSTF(1):ind.LSIG(2)),:); 
    lltfIdx = helperSymbolTiming(nonhtfields,conf.ChannelBandwidth);

    % Synchronize the received waveform given the offset between the
    % expected start of the L-LTF and actual start of L-LTF
    pktOffset = pktOffset+lltfIdx-double(ind.LLTF(1));
    % If no L-LTF detected or if packet detected outwith the range of
    % expected delays from the channel modeling; packet error
	if isempty(lltfIdx) || pktOffset<0
        rxPSDU = [];
        return;
    end        
    rxWaveform_awgn = rxWaveform_awgn(1+pktOffset:end,:);

    % Extract L-LTF and perform fine frequency offset correction
    lltf = rxWaveform_awgn(ind.LLTF(1):ind.LLTF(2),:); 
    fineFreqOff = wlanFineCFOEstimate(lltf,conf.ChannelBandwidth);
    rxWaveform_awgn = step(PFO,rxWaveform_awgn,-fineFreqOff);
    release(PFO); % Release object for subsequent processing

    % Estimate noise power in VHT fields
    lltf = rxWaveform_awgn(ind.LLTF(1):ind.LLTF(2),:); 
    demodLLTF = wlanLLTFDemodulate(lltf,conf.ChannelBandwidth);
    nVarVHT = helperNoiseEstimate(demodLLTF,conf);

    % Extract VHT-LTF samples from the waveform, demodulate and perform
    % channel estimation
    vhtltf = rxWaveform_awgn(ind.VHTLTF(1):ind.VHTLTF(2),:);
    vhtltfDemod = wlanVHTLTFDemodulate(vhtltf,conf);
    chanEst = wlanVHTLTFChannelEstimate(vhtltfDemod,conf);

    % Recover the transmitted PSDU in VHT Data
    % Extract VHT Data samples from the waveform and recover the PSDU
    vhtdata = rxWaveform_awgn(ind.VHTData(1):ind.VHTData(2),:);
    rxPSDU = wlanVHTDataRecover(vhtdata,chanEst,nVarVHT,conf);
end
end