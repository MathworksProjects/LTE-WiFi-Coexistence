function [BER, e_dct, pkt_ok, throughput, BER_new] = wiFi_report( conf, txPSDU, rxPSDU, off_wifi, off_wifi_dct, nWiFiFrames)
    BER = []; BER_new = [];
    e_dct = 0;
    pkt_ok = 0;

    tol = 10;               % Margin samples on the frame start detection
                            % due to multipath effect
    for fr = 1:nWiFiFrames
        %% Detection Check
        e_dct1 = 0;
        if(fr == 1)
            if(abs(off_wifi_dct(fr) - off_wifi) > tol);
                e_dct = e_dct + 1; 
                e_dct1 = 1;
            end
        else
            if(abs(off_wifi_dct(fr) - 0 ) > tol); 
                e_dct = e_dct + 1;
                e_dct1 = 1;
            end
        end

        %% BER Check
        e_ber = 0;
        for k = 1:conf.PSDULength*8
            if(txPSDU(k)~=rxPSDU(k,fr))
                e_ber = e_ber + 1;
            end
        end
        BER = [BER 100*e_ber/(conf.PSDULength*8)];              %#ok<AGROW>
        
        %% BER Conditional Check
        if(e_dct1==0)
            BER_new = [BER_new BER];                            %#ok<AGROW>
        else
            BER_new = [BER_new NaN];                            %#ok<AGROW>
        end

        %% Packet Success
        pkt_ok = pkt_ok + (1 - any(biterr(txPSDU,rxPSDU(:,fr)))) * (1 - e_dct1);
    end

    BER = mean(BER);
    
    %% Throughput
    throughput = (pkt_ok*length(txPSDU))/1e-3;
end