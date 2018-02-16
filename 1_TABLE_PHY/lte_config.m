%========================= IMPORTANT NOTE ================================%
% ** Author: Carlos Bocanegra Guerra, PhD student
% ** email: bocanegrac@coe.neu.edu
% ** Northeastern University, Boston MA
% 
% lte_config function configures the ABS Pattern, the channel and the
% channel estimatiton techniques carried out at the receiver.
%
% The ABS Pattern is configured through the parameter ABSConfig, which
% contains an array of the subframes that should not carry LTE data,
% allowin WiFi to schedule its transmissions during them. The user may add
% some configurations in this code as he needs it. Configurations 0 through
% 10 and 999 should not be erased from the system.
%
% The channel is defined by the number of taps, its average power and its
% delays. The channel estimation techniques are already defined by the
% Matlab LTE System Toolbox.
%
% [1] DownlinkChannelEstimationEqualizationExample.m
% 
%=========================================================================%

function [ ABSList, fsLTE, enbConfig, channel, cec] = lte_config(ABSConfig, nDLRB, cqiOffset, delayProfile)
    if    (ABSConfig == 0);   ABSList = 0;
    elseif(ABSConfig == 1);   ABSList = 1;
    elseif(ABSConfig == 2);   ABSList = 2;
    elseif(ABSConfig == 3);   ABSList = 3;
    elseif(ABSConfig == 4);   ABSList = 4;
    elseif(ABSConfig == 5);   ABSList = 5;
    elseif(ABSConfig == 6);   ABSList = 6;
    elseif(ABSConfig == 7);   ABSList = 7;
    elseif(ABSConfig == 8);   ABSList = 8;
    elseif(ABSConfig == 9);   ABSList = 9;
    elseif(ABSConfig == 10);  ABSList = [0,1];
	elseif(ABSConfig == 11);  ABSList = [1,2];
    elseif(ABSConfig == 999); ABSList = [];
    end

    if(any(cqiOffset==ABSList))
        fprintf('** ERROR: cqiOffset is included in ABSList');
        fprintf('Please chose a different value\n');
        return;
    end

    if(nDLRB == 6);    fsLTE = 1.92e6;           %specified in the standard
    elseif(nDLRB == 15);  fsLTE = 3.84e6;        %specified in the standard
    elseif(nDLRB == 25);  fsLTE = 7.68e6;        %specified in the standard
    elseif(nDLRB == 50);  fsLTE = 15.36e6;       %specified in the standard
    elseif(nDLRB == 75);  fsLTE = 23.04e6;       %specified in the standard
    elseif(nDLRB == 100);  fsLTE = 30.72e6;      %specified in the standard
    end

    enbConfig.NDLRB = nDLRB;                     % Number of resource blocks
    enbConfig.CellRefP = 1;                      % One transmit antenna port
    enbConfig.NCellID = 10;                      % Cell ID
    enbConfig.CyclicPrefix = 'Normal';           % Normal cyclic prefix
    enbConfig.DuplexMode = 'FDD';                % FDD
    enbConfig.PHICHDuration = 'Normal';          % Normal PHICH duration
    enbConfig.CFI = 3;                           % 4 PDCCH symbols
    enbConfig.Ng = 'Sixth';                      % HICH groups
    enbConfig.PDSCH.RVSeq = 0;                   % Disable HARQ
    enbConfig.PDSCH.CSIMode = 'PUCCH 1-0';       % Configure the CSI reporting mode
    enbConfig.PDSCH.Rho = 0;                     % DL Power allocation TS36.101 
                                                 % Table 8.2.1.1.1-1
    enbConfig.PDSCH.NLayers = 1;
    enbConfig.PDSCH.TxScheme = 'Port0';
    enbConfig.PDSCH.RNTI = 1;
    enbConfig.PDSCH.RV = 0;
    enbConfig.PDSCH.PRBSet = (0:enbConfig.NDLRB-1).';
    enbConfig.PDSCH.Modulation = 'QPSK';         % Initial Modulation

    %% CHANNEL CONFIGURATION
    %   MODEL A
    tauA = [0 10 20 30 40 50 60 70 80 90 110 140 170 200 240 290 340 390]*1e-9;
    pdbA = [0 -0.9 -1.7 -2.6 -3.5 -4.3 -5.2 -6.1 -6.9 -7.8 -4.7 -7.3 -9.9 -12.5 -13.7 -18.0 -22.4 -26.7];
    %   MODEL B
    tauB = [0 10 20 30 50 80 110 140 180 230 280 330 380 430 490 560 640 730]*1e-9;
    pdbB = [-2.6 -3 -3.5 -3.9 0 -1.3 -2.6 -3.9 -3.4 -5.6 -7.7 -9.9 -12.1 -14.3 -15.4 -18.4 -20.7 -24.6];
    %   MODEL C
    tauC = [0 10 20 30 50 80 110 140 180 230 280 330 400 490 600 730 880 1050]*1e-9;
    pdbC = [-3.3 -3.6 -3.9 -4.2 0 -0.9 -1.7 -2.6 -1.5 -3.0 -4.4 -5.9 -5.3 -7.9 -9.4 -13.2 -16.3 -21.2];
    %   MODEL D
    tauD = [0 10 20 30 50 80 110 140 180 230 280 330 400 490 600 730 880 1050]*1e-9;
    pdbD = [0 -10 -10.3 -10.6 -6.4 -7.2 -8.1 -9.0 -7.9 -9.4 -10.8 -12.3 -11.7 -14.3 -15.8 -19.6 -22.7 -27.6];
    %   MODEL E
    tauE = [0 10 20 40 70 100 140 190 240 320 430 560 710 880 1070 1280 1510 1760]*1e-9;
    pdbE = [-4.9 -5.1 -5.21122 -0.8 -1.3 -1.9 -0.3 -1.2 -2.1 -0.0 -1.9 -2.8 -5.4 -7.3 -10.6 -13.4 -17.4 -20.9];

    mapTau = containers.Map({'Model-A', 'Model-B', 'Model-C', 'Model-D', 'Model-E'},...
                            {tauA, tauB, tauC, tauD, tauE});
    mapPdb = containers.Map({'Model-A', 'Model-B', 'Model-C', 'Model-D', 'Model-E'},...
                            {pdbA, pdbB, pdbC, pdbD, pdbE});

    channel.Seed = randi(100);                          % Random channel seed
    channel.NRxAnts = 1;                                % 1 receive antenna
    channel.DelayProfile = 'Custom';                    % Delay spread (EPA, EVA, ETU, Custom)
    channel.MIMOCorrelation = 'Low';                    % Low (no) MIMO correlation
    channel.InitTime = 0;                               % Initialize at time zero
    channel.NTerms = 16;                                % Oscillators used in fading model
    channel.ModelType = 'GMEDS';                        % Rayleigh fading model type
    channel.InitPhase = 'Random';                       % Random initial phases
    channel.NormalizePathGains = 'On';                  % Normalize delay profile power 
    channel.NormalizeTxAnts = 'On';                     % Normalize for transmit antennas
    channel.DopplerFreq = 120;                          % Doppler frequency

    ofdmInfo = lteOFDMInfo(enbConfig);
    channel.SamplingRate = ofdmInfo.SamplingRate;       % Channel model sampling rate

    channel.PathDelays = mapTau(delayProfile);          % Channel Model TGn - Tap Delay
    channel.AveragePathGaindB = mapPdb(delayProfile);   % Channel Model TGn - Tap Gain

    %% CHANNEL ESTIMATOR CONFIGURATION
    cec.PilotAverage = 'UserDefined';                   % Pilot averaging method
    cec.FreqWindow = 9;                                 % Frequency averaging window in REs
    cec.TimeWindow = 9;                                 % Time averaging window in REs
    cec.InterpType = 'Cubic';                           % Cubic interpolation
    cec.InterpWinSize = 3;                              % Interpolate up to 3 subframes simultaneously
    cec.InterpWindow = 'Centred';                       % Interpolation windowing method
end

