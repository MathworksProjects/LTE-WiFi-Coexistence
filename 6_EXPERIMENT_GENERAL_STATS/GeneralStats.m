clear; close all; clc;

runnable = true;   %#ok % Always set to true if runnable

% Load the MAC Variables
run('../Configuration_MAC');

%% Global Experiment
exp = 'Nnodes'; varList = (1:1:20); alpha = 1;
% exp = 'Nnodes'; varList = [5 10 15 20]; alpha = 1;
% exp = 'Nnodes'; varList = [1 2 3 20]; alpha = 1;
% exp = 'Nnodes'; varList = [5 10 15 20]; alpha = 1;
% exp = 'PERth'; varList = (0:0.1:0.9); alpha = 1;
% exp = 'PSRth'; varList = (1:-0.05:0.1); alpha = 1;
% exp = 'BSloc'; varList = (10:5:120); alpha = 1;
% exp = 'Alpha'; varList = (0:0.1:1); alpha = 1;
GeneralStats_experiment(exp,varList,alpha, conf);

%% Closing Clause
runnable = false;

%% Function Global Experiment
function []= GeneralStats_experiment(experiment, varList, alpha, conf)
    Niter = conf.Niter;
    ABSList = [0 1];

	% Load MAC Tables (PSR) -> Average time to Transmit
    tableMACName = strcat(conf.TablePath,conf.TableName);
    load(tableMACName,conf.varMACList{:});
    conf.THPERLIST  = THPERLIST;
    conf.THPSDULIST = THPSDULIST;
    [~,idxAPtoGO]   = min(abs(repmat(alpha,size(ALPHAORG)) - ALPHAORG));
    [~,idxGOtoP2P]  = min(abs(repmat(1-alpha,size(ALPHAORG)) - ALPHAORG));
    conf.THTXTIMELIST_APtoGO  = THTXTIMELIST(idxAPtoGO,:);  %#ok % AP->GO
    conf.THTXTIMELIST_GOtoP2P = THTXTIMELIST(idxGOtoP2P,:);      % GO->P2P
    conf.THTXTIMELIST = THTXTIMELIST(end,:);                     % AP->P2P (alpha=1)

    % Load CONTENTION Tables -> Interference + Wi-Fi Contention!
    tableCONTENTIONName = strcat(conf.TablePath,conf.TableNameContention);
    load(tableCONTENTIONName,conf.varCONTENTList{:});
    conf.TXTIMECONTENT = TXTIMECONTENT;

    p1 = zeros(length(varList),4,length(ABSList));
    p2N = zeros(length(varList),2,length(ABSList));
    p4 = zeros(length(varList),4,length(ABSList));
    p5 = zeros(length(varList),4,length(ABSList));
    p6N = zeros(length(varList),4,length(ABSList));
    p7N = zeros(length(varList),4,length(ABSList));
    p8N = cell(5,length(varList),length(ABSList));
    p2 = []; p21 = []; p6 = []; p7 = []; p8 = []; p9 = cell(length(ABSList),1);
    p10 = zeros(length(varList),11,length(ABSList));
    for ABSCfg = ABSList
        y1 = []; y2 = []; y4 = []; y5 = []; y6 = []; y7 = []; y8 = [];
        y9 = cell(length(varList),1); y10 = []; y11 = [];

        % Load PSR/PER Tables (PRX,SINR) -> PSR
        if(ABSCfg == 0);     tablePSRName = strcat(conf.TablePath,'TABLE_PER_ABS0');
        elseif(ABSCfg == 1); tablePSRName = strcat(conf.TablePath,'TABLE_PER_ABS1');
        end
        load(tablePSRName,conf.varPSRList{:});
        conf.PERLIST  = PERLIST./100;
        conf.PRXLIST  = PRXLIST;
        conf.SINRLIST = SINRLIST;
        
    for var = varList
        if(strcmp(experiment,'Nnodes')); conf.Nnodes = var; fprintf('Nnodes = %d\n',conf.Nnodes);
        elseif(strcmp(experiment,'PERth')); conf.PER_th = var; conf.PSR_th = 1-conf.PER_th; fprintf('PER_th = %.2f\n',conf.PER_th);
        elseif(strcmp(experiment,'PSRth')); conf.PSR_th = var; fprintf('PSR_th = %.2f\n',conf.PSR_th);
        elseif(strcmp(experiment,'BSloc')); conf.posBS = [var 0]; fprintf('posBS = %d\n',conf.posBS(1));
        elseif(strcmp(experiment,'Alpha')); alpha = var; fprintf('Alpha = %.1f\n',alpha);
        end

        [~,idxAPtoGO]  = min(abs(repmat(alpha,size(ALPHAORG)) - ALPHAORG));
        [~,idxGOtoP2P] = min(abs(repmat(1-alpha,size(ALPHAORG)) - ALPHAORG));
        conf.THTXTIMELIST_APtoGO = THTXTIMELIST(idxAPtoGO,:);   % AP->GO 
        conf.THTXTIMELIST_GOtoP2P = THTXTIMELIST(idxGOtoP2P,:); % GO->P2P
        conf.THTXTIMELIST = THTXTIMELIST(end,:);                % AP->Rest

        rGO = zeros(Niter,1);
        rP2P = zeros(Niter,1);
        rClients_SZ = zeros(Niter,1);
        rClients_NSZ = zeros(Niter,1);
        ravPSR = zeros(Niter,1);
        ravPSR_new = zeros(Niter,1);
        ravThptGO = zeros(Niter,1);
        ravThptP2P = zeros(Niter,1);
        ravThptP2P_new = zeros(Niter,1);
        ravThptCSZ = zeros(Niter,1);
        ravThptCNSZ = zeros(Niter,1);
        K_avg = zeros(Niter,1);
        Ktot = [];
        totTxLTE_U1 = [];
        totTxEFi1 = [];
        parfor iter = 1:Niter
            [avThpt,avPSR,nodes,~,~,~,~,K_avg(iter),Ktot1,totTx] = EFi_function(conf);
            ravPSR(iter) = 100*avPSR.old;
            ravPSR_new(iter) = 100*avPSR.new;
            ravThptGO(iter) = avThpt.GO;
            ravThptP2P(iter) = avThpt.P2P;
            ravThptP2P_new(iter) = avThpt.P2P_new;
            ravThptCSZ(iter) = avThpt.CSZ;
            ravThptCNSZ(iter) = avThpt.CNSZ;
            Ktot = [Ktot Ktot1];    
            totTxLTE_U1 = [totTxLTE_U1 totTx.LTEU];
            totTxEFi1 = [totTxEFi1 totTx.EFi.'];

            rGO(iter) = 100 * length(nodes.GO) / conf.Nnodes;
            rP2P(iter) = 100 * length(nodes.P2P) / conf.Nnodes;
            rClients_SZ(iter) = 100 * length(nodes.Clients_SZ) / conf.Nnodes;
            rClients_NSZ(iter) = 100 * length(nodes.Clients_NSZ) / conf.Nnodes; 
        end
        ravPSR = ravPSR(~isnan(ravPSR));              % Get ride of NaN
        ravPSR_new = ravPSR_new(~isnan(ravPSR_new));  % Get ride of NaN

        y1 = [y1 ; mean(rGO) mean(rClients_SZ) mean(rClients_NSZ) mean(rP2P)];   %#ok<AGROW>
        y2 = [y2 ; mean(ravPSR) mean(ravPSR_new)];                               %#ok<AGROW>

        spIdx = find(ABSCfg==ABSList,1);
        p8N{1,var==varList,spIdx} = ravThptGO(~isnan(ravThptGO));         % Get ride of NaN
        p8N{2,var==varList,spIdx} = ravThptCSZ(~isnan(ravThptCSZ));               % Get ride of NaN
        p8N{3,var==varList,spIdx} = ravThptCNSZ(~isnan(ravThptCNSZ));            % Get ride of NaN
        p8N{4,var==varList,spIdx} = ravThptP2P(~isnan(ravThptP2P));               % Get ride of NaN
        p8N{5,var==varList,spIdx} = ravThptP2P_new(~isnan(ravThptP2P_new));   % Get ride of NaN

        y4 = [y4 ; mean(p8N{1,var==varList,spIdx}) mean(p8N{2,var==varList,spIdx}) mean(p8N{3,var==varList,spIdx}) mean(p8N{4,var==varList,spIdx})];      %#ok<AGROW>
        y5 = [y5 ; mean(p8N{1,var==varList,spIdx}) mean(p8N{2,var==varList,spIdx}) mean(p8N{3,var==varList,spIdx}) mean(p8N{5,var==varList,spIdx})];      %#ok<AGROW>

        ravThptGO2 = (mean(rGO)/100)*mean(p8N{1,var==varList,spIdx});
        ravThptCSZ2 = (mean(rClients_SZ)/100)*mean(p8N{2,var==varList,spIdx});
        ravThptCNSZ2 = (mean(rClients_NSZ)/100)*mean(p8N{3,var==varList,spIdx});
        ravThptP2P2 = (mean(rP2P)/100)*mean(p8N{4,var==varList,spIdx});
        ravThptP2P_new2 = (mean(rP2P)/100)*mean(p8N{5,var==varList,spIdx});

        y6 = [y6 ; ravThptGO2 ravThptCSZ2 ravThptCNSZ2 ravThptP2P2];      %#ok<AGROW>
        y7 = [y7 ; ravThptGO2 ravThptCSZ2 ravThptCNSZ2 ravThptP2P_new2];  %#ok<AGROW>
        
        y8 = [y8 mean(K_avg)];                                            %#ok<AGROW>
        y9{var==varList} = Ktot;
        
        y10 = [y10 mean(totTxLTE_U1)];                                    %#ok<AGROW>
        y11 = [y11 mean(totTxEFi1,2)];                                    %#ok<AGROW>

        fprintf('NGO      = %.2f%% ; NP2P     = %.2f%% ; NCSZ     = %.2f%% ; NCNSN     = %.2f%%\n',...
                mean(rGO), mean(rP2P), mean(rClients_SZ), mean(rClients_NSZ));
        fprintf('THGO     = %.2f(kbps) ; THP2P     = %.2f(kbps); THCSZ     = %.2f(kbps) ; THCNSN     = %.2f(kbps)\n',...
                mean(p8N{1,var==varList,spIdx}), mean(p8N{4,var==varList,spIdx}), mean(p8N{2,var==varList,spIdx}), mean(p8N{3,var==varList,spIdx}));
        fprintf('THGO     = %.2f(kbps) ; THP2P     = %.2f(kbps); THCSZ     = %.2f(kbps) ; THCNSN     = %.2f(kbps)\n',...
                mean(p8N{1,var==varList,spIdx}), mean(p8N{5,var==varList,spIdx}), mean(p8N{2,var==varList,spIdx}), mean(p8N{3,var==varList,spIdx}));
        fprintf('Avg. connections = %.3f\n',y8(end));
    %     fprintf('N=%d - avPSR=%.2f - avPSR_new=%.2f\n',Nnodes,mean(ravPSR),mean(ravPSR_new));
    end

    p2 = [p2 y2(:,2)];                      %#ok<AGROW>
    p21 = [p21 y2(1,1)];                    %#ok<AGROW>
	t = [y6(:,1) y6(:,2) y6(:,3) y6(:,4)];
    p6 = [p6 nansum(t.').'];                %#ok<AGROW>
    t = [y7(:,1) y7(:,2) y7(:,3) y7(:,4)];
    p7 = [p7 nansum(t.').'];                %#ok<AGROW>
    p8 = [p8 ; y8];                         %#ok<AGROW>
    p9{ABSCfg==ABSList} = y9;
    p10(:,:,ABSCfg==ABSList) = [y10;y11].';  % LTE-U and E-Fi

    % Data Node Categorization
    p1(:,:,find(ABSCfg==ABSList,1)) = y1;

    % Data Average PSR in % (802.11 V.S. E-Fi)
    p2N(:,:,find(ABSCfg==ABSList,1)) = y2;

    % Data the Average Throughput in Kbps for node categories (1)
    p4(:,:,find(ABSCfg==ABSList,1)) = y4;
    p5(:,:,find(ABSCfg==ABSList,1)) = y5;

    % Data the Average Throughput in Kbps for node categories (2)
    p6N(:,:,find(ABSCfg==ABSList,1)) = y6;
    p7N(:,:,find(ABSCfg==ABSList,1)) = y7;

    end

    % Save Data
    timeStamp = regexprep(string(datetime('now')),' ','-');  %#ok
    save('dataGeneralStats');

    % General plotting
    GeneralStatsPlot;
end