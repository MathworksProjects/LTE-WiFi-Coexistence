function [cqiBuffer, cqiBuffer_ABS, CQICalc, CQICalc_ABS, crc, offset, offsets] = ...
             lte_rx(sf, enbConfig, ABSList, rxWaveform_awgn, TBS, cqiBuffer, cqiBuffer_ABS, ...
             CQICalc, CQICalc_ABS, cqiCorrection, cqiPeriod, cqiOffset, cqiPtr, cec, offsets)
    %% LTE OFFSET CORRECTION
    if (mod(sf,10)==0)
        offset = lteDLFrameOffset(enbConfig,rxWaveform_awgn);  
        if (offset > 25)
            offset = offsets(end);
        end
    else
        offset = 0;
    end
    rxWaveform_awgn = rxWaveform_awgn(1+offset:end,:);

    %% LTE DEMODULATION
    rx_subframe = lteOFDMDemodulate(enbConfig,rxWaveform_awgn);

    %% CHANNEL ESTIMATION
    [chEstGrid,noiseEst] = lteDLChannelEstimate(enbConfig,enbConfig.PDSCH,cec,rx_subframe);

    %% CQI SELECTION
    [thisCQI,~] = lteCQISelect(enbConfig,enbConfig.PDSCH,chEstGrid,noiseEst);
    thisCQI = thisCQI - cqiCorrection;

    %% CQI UPDATE
    if( ~any(enbConfig.NSubframe==ABSList) )
        % CQI is UPDATED
        CQICalc = [CQICalc thisCQI];
        
        if(sf == cqiOffset)
            %CQI non-ABS subframes is UPDATED and REPORTED
            cqiBuffer(cqiPtr+1) = thisCQI;

            %CQI ABS subframes is NOT REPORTED
            cqiBuffer_ABS(cqiPtr+1) = cqiBuffer_ABS(mod(cqiPtr-1,cqiPeriod)+1);

        elseif(mod(sf,cqiPeriod)==cqiOffset)
            %CQI non-ABS subframes is REPORTED
            if(~isempty(CQICalc))
                cqiBuffer(cqiPtr+1) = ceil(sum(CQICalc)/length(CQICalc)); else
                cqiBuffer(cqiPtr+1) = cqiBuffer_ABS(mod(cqiPtr-1,cqiPeriod)+1); 
            end

            %CQI ABS subframes is REPORTED
            if(~isempty(CQICalc_ABS))
                cqiBuffer_ABS(cqiPtr+1) = ceil(sum(CQICalc_ABS)/length(CQICalc_ABS)); else
                cqiBuffer_ABS(cqiPtr+1) = cqiBuffer_ABS(mod(cqiPtr-1,cqiPeriod)+1);
            end
            
            CQICalc = [];
            CQICalc_ABS = [];
            
%             fprintf('-------------------- Update non-ABS CQI %d -> %d\n', CQI, cqiBuffer(cqiPtr+1));
%             fprintf('-------------------- Update     ABS CQI %d -> %d\n', CQI, cqiBuffer_ABS(cqiPtr+1));
        else
            % CQI non-ABS is NOT REPORTED
            cqiBuffer(cqiPtr+1) = cqiBuffer(mod(cqiPtr-1,cqiPeriod)+1);

            % CQI ABS is NOT REPORTED
            cqiBuffer_ABS(cqiPtr+1) = cqiBuffer_ABS(mod(cqiPtr-1,cqiPeriod)+1);
        end

        %% PDSCH DECODING
        ind = ltePDSCHIndices(enbConfig,enbConfig.PDSCH,enbConfig.PDSCH.PRBSet);                
        pdschRx = lteExtractResources(ind,rx_subframe) * (10^(-enbConfig.PDSCH.Rho/20));
        pdschChEst = lteExtractResources(ind,chEstGrid);
        [rxBits,~] = ltePDSCHDecode(enbConfig,enbConfig.PDSCH, ...
                                pdschRx,pdschChEst,noiseEst);

        % Decode the DL-SCH
        enbConfig.PDSCH.NTurboDecIts = 5;
        [~,crc] = lteDLSCHDecode(enbConfig,enbConfig.PDSCH,TBS,rxBits);     

        % Throughput calculation in subframe
        throughput = (TBS*(1-crc)/1e-2)*1e-6;

%         fprintf('nABS %d - %d - %3f \n',sf, CQI,(TBS*(1-crc)/1e-2)*1e-6);
    else
        % CQI is UPDATED
        CQICalc_ABS = [CQICalc_ABS thisCQI];

        if(mod(sf,cqiPeriod)==cqiOffset)
            %CQI non-ABS subframes is REPORTED
            if(~isempty(CQICalc))
                cqiBuffer(cqiPtr+1) = ceil(sum(CQICalc)/length(CQICalc)); else
                cqiBuffer(cqiPtr+1) = cqiBuffer_ABS(mod(cqiPtr-1,cqiPeriod)+1); 
            end

            %CQI ABS subframes is REPORTED
            if(~isempty(CQICalc_ABS))
                cqiBuffer_ABS(cqiPtr+1) = ceil(sum(CQICalc_ABS)/length(CQICalc_ABS)); else
                cqiBuffer_ABS(cqiPtr+1) = cqiBuffer_ABS(mod(cqiPtr-1,cqiPeriod)+1);
            end

            CQICalc = [];
            CQICalc_ABS = [];

%             fprintf('-------------------- Update non-ABS CQI %d -> %d\n', CQI, cqiBuffer(cqiPtr+1));
%             fprintf('-------------------- Update     ABS CQI %d -> %d\n', CQI, cqiBuffer_ABS(cqiPtr+1));
        else
            % CQI non-ABS is NOT UPDATED
            cqiBuffer(cqiPtr+1) = cqiBuffer(mod(cqiPtr-1,cqiPeriod)+1);

            % CQI ABS is NOT REPORTED
            cqiBuffer_ABS(cqiPtr+1) = cqiBuffer_ABS(mod(cqiPtr-1,cqiPeriod)+1);
        end
        
        % LTE Throughput is NOT UPDATED
        throughput = 0.0;
        
        % Subframe with no data modeled as an Error
        crc = 1;

%         fprintf('ABS  %d - %d - N/A \n',sf, CQI);
    end
end

