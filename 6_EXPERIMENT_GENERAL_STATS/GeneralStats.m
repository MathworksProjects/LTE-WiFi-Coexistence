clear; close all; clc;

runnable = true;   %#ok % Always set to true if runnable

% Load the MAC Variables
run('../Configuration_MAC');

%% Global Experiment
% exp = 'Nnodes'; varList = (1:1:20); alpha = 1;
% exp = 'Nnodes'; varList = [5 10 15 20]; alpha = 1;
% exp = 'Nnodes'; varList = [1 2 3 20]; alpha = 1;
exp = 'Nnodes'; varList = [5 10 15 20]; alpha = 1;
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
            Ktot = [Ktot Ktot1];                       %#ok<AGROW>
            totTxLTE_U1 = [totTxLTE_U1 totTx.LTEU];    %#ok<AGROW>
            totTxEFi1 = [totTxEFi1 totTx.EFi.'];       %#ok<AGROW>

            rGO(iter) = 100 * length(nodes.GO) / conf.Nnodes;
            rP2P(iter) = 100 * length(nodes.P2P) / conf.Nnodes;
            rClients_SZ(iter) = 100 * length(nodes.Clients_SZ) / conf.Nnodes;
            rClients_NSZ(iter) = 100 * length(nodes.Clients_NSZ) / conf.Nnodes; 
        end
        ravPSR = ravPSR(~isnan(ravPSR));              % Get ride of NaN
        ravPSR_new = ravPSR_new(~isnan(ravPSR_new));  % Get ride of NaN
        
        y1 = [y1 ; mean(rGO) mean(rClients_SZ) mean(rClients_NSZ) mean(rP2P)];   %#ok<AGROW>
        y2 = [y2 ; mean(ravPSR) mean(ravPSR_new)];                               %#ok<AGROW>

        ravThptGO1 = ravThptGO(~isnan(ravThptGO));                  % Get ride of NaN
        ravThptP2P1 = ravThptP2P(~isnan(ravThptP2P));               % Get ride of NaN
        ravThptP2P_new1 = ravThptP2P_new(~isnan(ravThptP2P_new));   % Get ride of NaN
        ravThptCSZ1 = ravThptCSZ(~isnan(ravThptCSZ));               % Get ride of NaN
        ravThptCNSZ1 = ravThptCNSZ(~isnan(ravThptCNSZ));            % Get ride of NaN
        
        y4 = [y4 ; mean(ravThptGO1) mean(ravThptCSZ1) mean(ravThptCNSZ1) mean(ravThptP2P1)];      %#ok<AGROW>
        y5 = [y5 ; mean(ravThptGO1) mean(ravThptCSZ1) mean(ravThptCNSZ1) mean(ravThptP2P_new1)];  %#ok<AGROW>

        ravThptGO2 = (mean(rGO)/100)*mean(ravThptGO1);
        ravThptCSZ2 = (mean(rClients_SZ)/100)*mean(ravThptCSZ1);
        ravThptCNSZ2 = (mean(rClients_NSZ)/100)*mean(ravThptCNSZ1);
        ravThptP2P2 = (mean(rP2P)/100)*mean(ravThptP2P1);
        ravThptP2P_new2 = (mean(rP2P)/100)*mean(ravThptP2P_new1);

        y6 = [y6 ; ravThptGO2 ravThptCSZ2 ravThptCNSZ2 ravThptP2P2];      %#ok<AGROW>
        y7 = [y7 ; ravThptGO2 ravThptCSZ2 ravThptCNSZ2 ravThptP2P_new2];  %#ok<AGROW>
        
        y8 = [y8 mean(K_avg)];                                            %#ok<AGROW>
        y9{var==varList} = Ktot;
        
        y10 = [y10 mean(totTxLTE_U1)];                                    %#ok<AGROW>
        y11 = [y11 mean(totTxEFi1,2)];                                    %#ok<AGROW>

        fprintf('NGO      = %.2f%% ; NP2P     = %.2f%% ; NCSZ     = %.2f%% ; NCNSN     = %.2f%%\n',...
                mean(rGO), mean(rP2P), mean(rClients_SZ), mean(rClients_NSZ));
        fprintf('THGO     = %.2f(kbps) ; THP2P     = %.2f(kbps); THCSZ     = %.2f(kbps) ; THCNSN     = %.2f(kbps)\n',...
                mean(ravThptGO1), mean(ravThptP2P1), mean(ravThptCSZ1), mean(ravThptCNSZ1));
        fprintf('THGO     = %.2f(kbps) ; THP2P     = %.2f(kbps); THCSZ     = %.2f(kbps) ; THCNSN     = %.2f(kbps)\n',...
                mean(ravThptGO1), mean(ravThptP2P_new1), mean(ravThptCSZ1), mean(ravThptCNSZ1));
        fprintf('Avg. connections = %.3f\n',y8(end));
    %     fprintf('N=%d - avPSR=%.2f - avPSR_new=%.2f\n',Nnodes,mean(ravPSR),mean(ravPSR_new));
    end

    labels = varList;
    if(strcmp(experiment,'Nnodes')); xslabel = 'Number of Nodes in the Area';
    elseif(strcmp(experiment,'PERth')); xslabel = 'PER Threshold for Safe Zone (%)'; labels = varList.*100;
	elseif(strcmp(experiment,'PSRth')); xslabel = 'PSR Threshold for Safe Zone (%)'; labels = varList.*100;
    elseif(strcmp(experiment,'BSloc')); xslabel = 'AP-BS distance (m)';
    elseif(strcmp(experiment,'Alpha')); xslabel = 'Alpha';
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
    
    % Plot Statistical Node Categorization
    figure(1);
    if(ABSCfg == 0)
        subplot(2,1,1)
    else
        subplot(2,1,2)
    end
    b = bar(y1,'stacked');
    b(4).FaceColor = [192 192 192]./255;
    b(3).FaceColor = [144 144 144]./255;
    b(2).FaceColor = [104 104 104]./255;
    b(1).FaceColor = [48 48 48]./255;
    xlim([0.5 length(labels)+0.5]); ylim([0 100]);
    set(gca, 'XTick', 1:length(labels)); % Change x-axis ticks
    set(gca, 'XTickLabel', labels); % Change x-axis ticks labels.
%     xlim([min(labels)-0.5 length(labels)+0.5]); ylim([0 100]);
    hxlab = xlabel(xslabel);
    hylab = ylabel('Nodes (%)');
%     hleg1  = legend(b(1:2),'Group Owners','S.Z. Clients');
%     ah=axes('position',get(gca,'position'),'visible','off');
%     hleg2  = legend(ah,b(3:4),'N.S.Z. Clients','Wi-Fi Direct Clients');
    hleg  = legend('Group Owners','S.Z. Clients','N.S.Z. Clients','Wi-Fi Direct Clients');
    str = sprintf('(ABS %d)',ABSCfg);
%     str = '';
    tit = title(str);
    set(tit,'FontSize',12);
    set(hxlab,'FontSize',12);
    set(hylab,'FontSize',12);
%     set(hleg,'FontSize',12,'Location','Northwest');
    set(hleg,'FontSize',12,'Location','Northwest');
    grid minor;

    % Plot Average PSR in % (802.11 V.S. E-Fi)
    figure(2)
    if(ABSCfg == 0)
        subplot(1,2,1)
    else
        subplot(1,2,2)
    end
%     b = bar(y2);
%     b(1).FaceColor = [104 104 104]./255;
%     b(2).FaceColor = [48 48 48]./255;
%     set(gca, 'XTick', 1:length(labels)); % Change x-axis ticks
%     set(gca, 'XTickLabel', labels); % Change x-axis ticks labels.
%     xlim([0 length(labels)+1]);
%     ylim([min(y2(:,1)-5) max(y2(:,2)+5)]);
    plot(varList,y2(:,1),'k-.','Marker','x','LineWidth',2);
    hold on;
    plot(varList,y2(:,2),'k-','Marker','s','LineWidth',2);
    xlim([min(varList) max(varList)]);
    ylim([min(y2(:,1)-5) max(y2(:,2)+5)]);
    hxlab = xlabel(xslabel);
    hylab = ylabel('Average PSR (%)');
    hleg  = legend('No Optimization','Modofied Hungarian');
    str = sprintf('(ABS %d) Average PSR',ABSCfg);
    tit = title(str);
    set(tit,'FontSize',12);
    set(hxlab,'FontSize',12);
    set(hylab,'FontSize',12);
    set(hleg,'FontSize',12);
    set(hleg,'Location','Northwest');
    grid minor;

    % Plot the Average Throughput in Kbps for node categories (1)
    figID = 3; figure(figID);
    if(ABSCfg == 0)
        subplot(2,1,1)
    else
        subplot(2,1,2)
    end
    NumGroupsPerAxis = length(labels);   % Number bars for x axis
    NumStacksPerGroup = 2;               % Comparing Original with Modified Hungarian
    NumStackElements = 4;                % 4 groups: GO, P2P, CSZ and CNSZ
    stackData = zeros(NumGroupsPerAxis,NumStacksPerGroup,NumStackElements);
    stackData(:,1,:) = y4;
    stackData(:,2,:) = y5;
    groupLabels = {labels};
    h = plotBarStackGroups(stackData, groupLabels,figID);
    hxlab = xlabel(xslabel);
    hylab = ylabel('Average Throughput (Kbps)');
    hleg1  = legend(h(1,:),'Group Owners','S.Z. Clients','N.S.Z. Clients','Wi-Fi Direct Clients');
    ah=axes('position',get(gca,'position'),'visible','off');
    hleg2  = legend(ah,h(2,:),'Group Owners','S.Z. Clients','N.S.Z. Clients','Wi-Fi Direct Clients');
    str = sprintf('(ABS %d) Av. User thpt. and contribution by categories to total Av. Thpt.',ABSCfg);
    tit = title(str);
    set(tit,'FontSize',12);
    set(hxlab,'FontSize',12);
    set(hylab,'FontSize',12);
    set(hleg1,'FontSize',12,'Location','Northwest');
    set(hleg2,'FontSize',12,'Location','Northeast');
    grid minor;

    % Plot the Average Throughput in Kbps for node categories (2)
    figID = 4; figure(figID);
    if(ABSCfg == 0)
        subplot(2,1,1)
    else
        subplot(2,1,2)
    end
    NumGroupsPerAxis = length(labels);   % Number bars for x axis
    NumStacksPerGroup = 2;               % Comparing Original with Modified Hungarian
    NumStackElements = 4;                % 4 groups: GO, P2P, CSZ and CNSZ
    stackData = zeros(NumGroupsPerAxis,NumStacksPerGroup,NumStackElements);
    stackData(:,1,:) = y6;
    stackData(:,2,:) = y7;
    groupLabels = {labels};
    h = plotBarStackGroups(stackData, groupLabels,figID);
    hxlab = xlabel(xslabel);
    hylab = ylabel('Average Throughput (Kbps)');
	hleg1  = legend(h(1,:),'Group Owners','S.Z. Clients','N.S.Z. Clients','Wi-Fi Direct Clients');
    ah=axes('position',get(gca,'position'),'visible','off');
    hleg2  = legend(ah,h(2,:),'Group Owners','S.Z. Clients','N.S.Z. Clients','Wi-Fi Direct Clients');
    str = sprintf('(ABS %d) Av. User thpt. and contribution by categories to total Av. Thpt.',ABSCfg);
    tit = title(str);
    set(hxlab,'FontSize',12);
    set(hylab,'FontSize',12);
    set(tit,'FontSize',12);
	set(hleg1,'FontSize',12,'Location','Northwest');
    set(hleg2,'FontSize',12,'Location','Northeast');
    grid minor;
    
    % Plot CDF of Throughput for node categories
    figID = 5; figure(figID);
    if(ABSCfg == 0)
        subplot(2,1,1)
    else
        subplot(2,1,2)
    end
    hold on
    h = cdfplot(ravThptGO1./1e3);
    %     h.Color = [53, 42, 134] / 255;
        h.Color = [0 0 0]./255;
        h.LineStyle = '--';
        h.LineWidth = 3;
    h = cdfplot(ravThptCSZ1./1e3);
    %     h.Color = [6, 156, 207] / 255;
        h.Color = [105 105 105]./255;
        h.LineStyle = '--';
        h.LineWidth = 3;
    if ~isempty(ravThptCNSZ1)
        h = cdfplot(ravThptCNSZ1./1e3);
        %     h.Color = [165, 190, 106] / 255;
        h.Color = [192 192 192]./255;
        h.LineStyle = ':';
        h.LineWidth = 3;
    end
    h = cdfplot(ravThptP2P1./1e3);
    %     h.Color = [248, 250, 3] / 255;
        h.Color = [144 144 144]./255;
        h.LineStyle = '-.';
        h.LineWidth = 3;
    h = cdfplot(ravThptP2P_new1./1e3);
    %     h.Color = [248, 134, 3] / 255;
        h.Color = [144 144 144]./255;
        h.LineWidth = 3;
	if ~isempty(ravThptCNSZ1)
        hl = legend('Group Owner', 'S.Z. Client', 'N.S.Z. Client',  'Wi-Fi Direct Client (Wi-Fi)', 'Wi-Fi Direct Client (E-Fi)');
    else
        hl = legend('Group Owner', 'S.Z. Client', 'Wi-Fi Direct Client (Wi-Fi)', 'Wi-Fi Direct Client (E-Fi)');
    end
    hl.FontSize = 12;
    hl.Location = 'southeast';
    hx = xlabel('Throughput (Mbps)');
    hx.FontSize = 12;
    hy = ylabel('CDF');
    hy.FontSize = 12;
    set(gca,'fontsize',12)
    str = sprintf('(ABS %d) Throughput CDF',ABSCfg);
%     str = '';
    ht = title(str);
    ht.FontSize = 12;
    hold off
    grid minor
    
    end

    % Plot Final Average PSR distribution in %
    figure(6)
    b = bar(p2);
    b(1).FaceColor = [104 104 104]./255;
    b(2).FaceColor = [48 48 48]./255;
    hold on
    tt(1) = plot((1:length(labels)),p21(1,1)*ones(length(labels),1),'--o','Color','k','LineWidth',1.5);
    tt(2) = plot((1:length(labels)),p21(1,2)*ones(length(labels),1),'--x','Color','k','LineWidth',1.5);
    set(gca, 'XTick', 1:length(labels)); % Change x-axis ticks
    set(gca, 'XTickLabel', labels); % Change x-axis ticks labels.
    grid minor;
    xlim([0.5 length(labels)+0.5]);
    ylim([min(p21(:,1))-5 max(p2(:,2)+8)]);
    ylim([min(p21(:,1))-5 max(p2(:,2)+20)]);
    ylim([min(p21(:,1))-5 max(p2(:,2)+8)]);
    hxlab = xlabel(xslabel);
    hylab = ylabel('Average PSR (%)');
    hleg1  = legend(b,'ABS0 Modified Hungarian','ABS1 Modified Hungarian');
    ah=axes('position',get(gca,'position'),'visible','off');
    hleg2  = legend(ah,tt,'ABS0 No Optimization','ABS1 No Optimization');
%     str = sprintf('Average PSR');
    str = '';
    tit = title(str);
    set(tit,'FontSize',14);
    set(hxlab,'FontSize',14);
    set(hylab,'FontSize',14);
    set(hleg1,'FontSize',12,'Location','Northwest');
    set(hleg2,'FontSize',12,'Location','Northeast');
    
    % Plot Final Average Throughput in Kbps
    figure(7)
    b = bar(p7);
    b(1).FaceColor = [104 104 104]./255;
    b(2).FaceColor = [48 48 48]./255;
    hold on
    plot((1:length(labels)),mean(p6(:,1))*ones(length(labels),1),'-o','Color','k')
    plot((1:length(labels)),mean(p6(:,2))*ones(length(labels),1),'-*','Color','k')
    set(gca, 'XTick', 1:length(labels)); % Change x-axis ticks
    set(gca, 'XTickLabel', labels); % Change x-axis ticks labels.
    xlim([0.5 length(labels)+0.5]);
    ylim([min(p6(:,1))-200 max(p6(:,2)+500)]);
    hxlab = xlabel(xslabel);
    hylab = ylabel('Average Throughput (Kbps)');
    hleg  = legend('ABS0 Modified Hungarian','ABS1 Modified Hungarian','ABS0 No Optimization','ABS1 No Optimization');
%     str = sprintf('Av. User thpt. and contribution by categories to total Av. Thpt.');
    str = '';
    tit = title(str);
    set(tit,'FontSize',14);
    set(hxlab,'FontSize',14);
    set(hylab,'FontSize',14);
    set(hleg,'FontSize',14);
    set(hleg,'Location','Northwest');
    grid minor;

    % Plot Final average number of WiFi clients per cluster
    figure(8); hold on;
    plot(varList,p8(1,:),'k--','LineWidth',2);
    plot(varList,p8(2,:),'k-','LineWidth',2);
    xlabel('Number of Wi-Fi nodes','FontSize',14);
    ylabel('Nodes','FontSize',14);
    title('Average nodes per Wi-Fi Direct group','FontSize',14)
    lg = legend('ABS 0','ABS 1');
    set(lg,'FontSize',14);
    grid minor
    
    % ECDF of the number of WiFi clients per cluster
    strleg = cell(length(varList),1);
    figure(9); hold on;
    tp = (ABSList==1);
    t = p9{tp};
    totT = [];
    for n = 1:length(varList)
        H = histcounts(t{n});
        H = 100.*H./sum(H);
        Htrim = H(1:5);
        totT = [totT ; Htrim];  %#ok
        strleg{n} = strcat('Nnodes=',num2str(varList(n)));
    end
    b = bar(totT.');
	b(4).FaceColor = [192 192 192]./255;
    b(3).FaceColor = [144 144 144]./255;
    b(2).FaceColor = [104 104 104]./255;
    b(1).FaceColor = [48 48 48]./255;
    xlabel('Number of Wi-Fi nodes','FontSize',14);
    ylabel('Percentage (%)','FontSize',14);
    title('Number of nodes per group - ABS1','FontSize',14);
    lg = legend(strleg);
    set(lg,'FontSize',14);
    grid minor
    
    % Latency Analysis . A comparison between LTE-U and E-Fi
    % Only plot Results for Configuration ABSCfg=1 (ABS 1)
    figure(10);
    b = bar(p10(:,:,ABSList==1));
    grid minor
    xlabel('Number of Nodes in the network','FontSize',14);
    ylabel('Required time (ms)','FontSize',14);
    b(1).FaceColor = 'r';
    hleg1 = legend(b(1),'LTE-U');
    ah=axes('position',get(gca,'position'),'visible','off');
    hleg2 = legend(ah,b(2:6),'1 ABS','2 ABS','3 ABS','4 ABS','5 ABS');
    ah=axes('position',get(gca,'position'),'visible','off');
    hleg3 = legend(ah,b(7:end),'6 ABS','7 ABS','8 ABS','9 ABS','10 ABS');
    set(hleg1,'FontSize',12,'Location','NorthWest');
    set(hleg2,'FontSize',12,'Location','SouthEast');
    set(hleg3,'FontSize',12,'Location','NorthEast');
    for i = 2:11
    b(i).FaceColor = [26*(i-2) 26*(i-2) 26*(i-2)]./255;
    end
    set(gca,'XTickLabel',[5 10 15 20] );
    title('Latency comparison - LTE-U vs E-Fi','FontSize',15);
    title('LTE-U vs E-Fi - Latency Analysis','FontSize',15);
end