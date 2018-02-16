function [txWaveform, subframe, enbConfig, TBS, CQI, cqiPtr] = lte_tx(sf, enbConfig, ABSList, cqiPeriod, cqiBuffer, Ptx_LTE)
    % Set subframe and frame number
    enbConfig.NSubframe = mod(sf,10);
    enbConfig.NFrame = floor(sf/10);

    % Generate empty subframe
    subframe = lteDLResourceGrid(enbConfig);

    cqiPtr = mod(sf,cqiPeriod);
    
    %% PDSCH
    if( ~any(enbConfig.NSubframe==ABSList) )
        % Select CQI, reading the oldest value from the CQI buffer        
        CQI = cqiBuffer(cqiPtr+1);

        % Select MCS according to CQI using TS36.101 Table A.4-1 CSI RMC
        % RC.1 FDD (MCS.1), which defines the relationship between CQI
        % indices and MCS indices
        IMCSTable = [-1 0 0 2 4 6 8 11 13 16 18 21 23 25 27 27];
        IMCS = IMCSTable(CQI+1);

        % Determine TBS and modulation order
        [ITBS,modulation] = lteMCS(IMCS);
        enbConfig.PDSCH.Modulation = {modulation};

        TBS = double(lteTBS(size(enbConfig.PDSCH.PRBSet,1),ITBS));
        enbConfig.PDSCH.TrBlkSizes(enbConfig.NSubframe+1) = TBS;

        [pdschIndices,pdschInfo] = ltePDSCHIndices(enbConfig, enbConfig.PDSCH, enbConfig.PDSCH.PRBSet, {'1based'});

        codedTrBlkSize = pdschInfo.G;
        dlschTransportBlk = randi([0 1], TBS, 1);
        codedTrBlock = lteDLSCH(enbConfig, enbConfig.PDSCH, codedTrBlkSize, dlschTransportBlk);

        pdschSymbols = ltePDSCH(enbConfig, enbConfig.PDSCH, codedTrBlock);

        subframe(pdschIndices) = pdschSymbols; 
    else
        TBS = 0; CQI = [];
    end
    
    %% PDCCH
    % PDCCH - DCI message configuration
    dciConfig.DCIFormat = 'Format1A';     % DCI message format
    dciConfig.Allocation.RIV = 26;        % Resource indication value

    [~, dciMessageBits] = lteDCI(enbConfig, dciConfig); % DCI message

    % PDCCH - DCI Channel Coding
    pdcchConfig.NDLRB = enbConfig.NDLRB;  % Number of DL-RB in total BW
    pdcchConfig.RNTI = 1;                 % 16-bit value number
    pdcchConfig.PDCCHFormat = 0;          % 1-CCE of aggregation level 1

    % PDCCH - DCI message bits coding
    codedDciBits = lteDCIEncode(pdcchConfig, dciMessageBits);

    % PDCCH - Bits Generation
    pdcchInfo = ltePDCCHInfo(enbConfig);    % Get the total resources for PDCCH
    pdcchBits = -1*ones(pdcchInfo.MTot, 1); % Initialized with -1

    candidates = ltePDCCHSpace(enbConfig, pdcchConfig, {'bits','1based'});
    pdcchBits( candidates(1, 1) : candidates(1, 2) ) = codedDciBits;

    pdcchSymbols = ltePDCCH(enbConfig, pdcchBits);
    pdcchIndices = ltePDCCHIndices(enbConfig, {'1based'});

    subframe(pdcchIndices) = pdcchSymbols;
    
    %% PSCH & SSCH 
    pssSymbols = ltePSS(enbConfig);
    sssSymbols = lteSSS(enbConfig);
    pssIndices = ltePSSIndices(enbConfig);
    sssIndices = lteSSSIndices(enbConfig);

    subframe(pssIndices) = pssSymbols;
    subframe(sssIndices) = sssSymbols;

    %% RS
    cellRsSym = lteCellRS(enbConfig);
    cellRsInd = lteCellRSIndices(enbConfig);

    subframe(cellRsInd) = cellRsSym;

    %% PCFICH
    cfiBits = lteCFI(enbConfig);
    pcfichSymbols = ltePCFICH(enbConfig, cfiBits);
    pcfichIndices = ltePCFICHIndices(enbConfig);

    subframe(pcfichIndices) = pcfichSymbols;

    %% PHICH
    phichSymbols = ltePHICH(enbConfig,[0,2,0]);
    phichIndices = ltePHICHIndices(enbConfig);
    
    subframe(phichIndices) = phichSymbols;
    
    %% PBCCH
    mib = lteMIB(enbConfig);
    bchCoded = lteBCH(enbConfig,mib);
    
    pbchIndices = ltePBCHIndices(enbConfig);
    pbchSymbols = ltePBCH(enbConfig,bchCoded);
    
    % The PBCH is always divided into 4 subframes (4 blocks of 240 complex 
    % symbols). It takes 4 frames to transmit the MIB, the PBCH is always
    % transmitted in the subframe 0 of every frame.
    if(~isempty(pbchIndices))
        subframe(pbchIndices) = ...
         pbchSymbols(240*mod(floor(sf/10),4)+1:240*(mod(floor(sf/10),4)+1));
    end

    %% OFDM MODULATION
    [txWaveform,~] = lteOFDMModulate(enbConfig,subframe);
    
    %% Transmitted Power
	Ptx_LTE = 10^((Ptx_LTE-30)/10);
    sym = 4;    % Reference Signals are sent and can be used to calculate
                % the Received LTE Power. They are located in symbols
                % 0, 4, 7 and 11
    lLTE = length(txWaveform);
    lte_ini_power = 1+sym*floor(lLTE/14);
    lte_end_power = (sym+1)*floor(lLTE/14);
    P_LTE  = coex_power(txWaveform(lte_ini_power:lte_end_power,1));
    txWaveform = sqrt(Ptx_LTE/P_LTE).*txWaveform;
end

