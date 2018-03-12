conf.Radius = 60;             % Coverage area where devices are deployed
conf.Cover_angle = pi/2;      % Either pi/2 (Half Area) or pi (Whole area)
conf.posAP = [0 0];           % AP Location in the XY plane
conf.posBS = [conf.Radius 0]; % BS Location in the XY plane
conf.Ptx_WiFi = 17;           % WiFi-AP transmit power
conf.Ptx_LTE = 17;            % LTE-BS transmit power
conf.Pnoise = -95;            % Noise Power in dBm
conf.Nnodes = 20;             % Number of WiFi Devices
conf.PER_th = 0.6;            % PER threshold defining the Safe Zone (0.4 CDF Th)
conf.PSR_th = 1-conf.PER_th;  % PSR threshold defining the Safe Zone
conf.K = 7;                   % Maximum Number of WiFi-Direct Connections
conf.PSDULength = 600;        % Number of bits in the Payload
conf.ABSCfg = 1;              % Possible Values: 0, 1, 5
conf.Criteria = 'Th_or';      % Possible Values: 'PER_or' or 'Th_or'
conf.MCS = 0;                 % MCS defines the max bits per packet
conf.Alpha = 1;               % Modified Congestion Window (Default is 1)
conf.PSR_minth = 0.9*conf.PSR_th;  % Tolerance for Feas. Region Exp. (Lower)
conf.PSR_maxth = 1.1*conf.PSR_th;  % Tolerance for Feas. Region Exp. (Upper)
conf.NS3Bridge = false;       % Generates the Bridge for NS-3
conf.NS3FileName = 'NS3Input_K4.txt'; % Name of the NS-3 Bridge file

conf.varPSRList = {'PERLIST','PRXLIST','SINRLIST'};
conf.varMACList = {'THTXTIMELIST','THPERLIST','THPSDULIST','ALPHAORG'};
conf.varCONTENTList = {'TXTIMECONTENT'};

conf.TablePath = '../0_DATA/';
% conf.TableName = 'TABLE_NEW_MAC';
conf.TableName = 'TABLE_NEW_MAC_testing';
conf.TableNameContention = 'TABLE_AV_TX_CONTENTION';

conf.blackAndWhite = false;

% Number of Iterations for the set of experiments. For stable and
% meaningful results, the following values are recommended:
% - CreateMacTable:            N/A
% - Categorization_experiment: 200 iterations
% - Throughput_experiment:     10000 iterations
% - FeasibleRegion_experiment: N/A
% - General_stats:             1000 iterations
conf.Niter = 1000;