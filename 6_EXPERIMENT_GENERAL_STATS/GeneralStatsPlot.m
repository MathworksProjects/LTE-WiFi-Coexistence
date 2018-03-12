clear; close all; clc;

% CONFIGURATION PARAMETERS %
conf.blackAndWhite = false;  % Plots in Grey Scale
% ------------------------ %

% Load variables from previous execution
load('dataGeneralStats');
fprintf('Plotting execution from %s\n',timeStamp);

labels = varList;
if(strcmp(experiment,'Nnodes')); xslabel = 'Number of Nodes in the Area';
elseif(strcmp(experiment,'PERth')); xslabel = 'PER Threshold for Safe Zone (%)'; labels = varList.*100;
elseif(strcmp(experiment,'PSRth')); xslabel = 'PSR Threshold for Safe Zone (%)'; labels = varList.*100;
elseif(strcmp(experiment,'BSloc')); xslabel = 'AP-BS distance (m)';
elseif(strcmp(experiment,'Alpha')); xslabel = 'Alpha';
end

for ABSCfg = ABSList
    % Get subplot index
    spIdx = find(ABSCfg==ABSList,1);
    
    % Plot Statistical Node Categorization
    figure(1); 
    subplot(2,1,spIdx);
    b = bar(p1(:,:,spIdx),'stacked');
    if conf.blackAndWhite
        b(4).FaceColor = [192 192 192]./255;
        b(3).FaceColor = [144 144 144]./255;
        b(2).FaceColor = [104 104 104]./255;
        b(1).FaceColor = [48 48 48]./255;
    else
        b(4).FaceColor = 'r';
        b(3).FaceColor = 'g';
        b(2).FaceColor = 'b';
        b(1).FaceColor = 'k';
    end
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
    figure(2);
    subplot(2,1,spIdx);
    % b = bar(p2N(:,:,spIdx));
    % b(1).FaceColor = [104 104 104]./255;
    % b(2).FaceColor = [48 48 48]./255;
    % set(gca, 'XTick', 1:length(labels)); % Change x-axis ticks
    % set(gca, 'XTickLabel', labels); % Change x-axis ticks labels.
    % xlim([0 length(labels)+1]);
    % ylim([min(p2N(:,1,spIdx)-5) max(p2N(:,2,spIdx)+5)]);
    plot(varList,p2N(:,1,spIdx),'k-.','Marker','x','LineWidth',2);
    hold on;
    plot(varList,p2N(:,2,spIdx),'k-','Marker','s','LineWidth',2);
    xlim([min(varList) max(varList)]);
    ylim([min(p2N(:,1,spIdx)-5) max(p2N(:,2,spIdx)+5)]);
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
    subplot(2,1,spIdx);
    NumGroupsPerAxis = length(labels);   % Number bars for x axis
    NumStacksPerGroup = 2;               % Comparing Original with Modified Hungarian
    NumStackElements = 4;                % 4 groups: GO, P2P, CSZ and CNSZ
    stackData = zeros(NumGroupsPerAxis,NumStacksPerGroup,NumStackElements);
    stackData(:,1,:) = p4(:,:,spIdx);
    stackData(:,2,:) = p5(:,:,spIdx);
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
    subplot(2,1,spIdx);
    NumGroupsPerAxis = length(labels);   % Number bars for x axis
    NumStacksPerGroup = 2;               % Comparing Original with Modified Hungarian
    NumStackElements = 4;                % 4 groups: GO, P2P, CSZ and CNSZ
    stackData = zeros(NumGroupsPerAxis,NumStacksPerGroup,NumStackElements);
    stackData(:,1,:) = p6N(:,:,spIdx);
    stackData(:,2,:) = p7N(:,:,spIdx);
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
    subplot(2,1,spIdx);
    hold on
    h1 = cdfplot(p8N{1,end,spIdx}./1e3);
    h1.LineWidth = 3;
    h2 = cdfplot(p8N{2,end,spIdx}./1e3);
    h2.LineWidth = 3;
    if ~isempty(p8N{3,end,spIdx})
        h3 = cdfplot(p8N{3,end,spIdx}./1e3);
        h3.LineWidth = 3;
    end
    h4 = cdfplot(p8N{4,end,spIdx}./1e3);
    h4.LineWidth = 3;
    h5 = cdfplot(p8N{5,end,spIdx}./1e3);
    h5.LineWidth = 3;
        
    if conf.blackAndWhite
        h1.Color = [0 0 0]./255;
        h2.Color = [105 105 105]./255;
        h3.Color = [192 192 192]./255;
        h4.Color = [144 144 144]./255;
        h5.Color = [144 144 144]./255;
        
        h1.LineStyle = '--';
        h2.LineStyle = '--';
        h3.LineStyle = ':';
        h4.LineStyle = '-.';
    else
        h1.Color = [53, 42, 134] / 255;
        h2.Color = [6, 156, 207] / 255;
        h3.Color = [165, 190, 106] / 255;
        h4.Color = [248, 250, 3] / 255;
        h5.Color = [248, 134, 3] / 255;
        
        h1.Color = 'k';
        h2.Color = 'b';
        h3.Color = 'r';
        h4.Color = 'g';
        h4.LineStyle = '--';
        h5.Color = 'g';
    end

	if ~isempty(p8N{3,end,spIdx})
        hl = legend('Group Owner', 'S.Z. Client', 'N.S.Z. Client',  'Wi-Fi Direct Client (Wi-Fi)', 'Wi-Fi Direct Client (E-Fi)');
    else
        hl = legend('Group Owner', 'S.Z. Client', 'Wi-Fi Direct Client (Wi-Fi)', 'Wi-Fi Direct Client (E-Fi)');
    end
    hl.FontSize = 12;
    set(hl,'FontSize',10,'Location','Southeast');
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
if conf.blackAndWhite
    b(1).FaceColor = [104 104 104]./255;
    b(2).FaceColor = [48 48 48]./255;
else
    b(1).FaceColor = [0 128 255]./255;
    b(2).FaceColor = [255 128 0]./255;
end
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
if conf.blackAndWhite
    b(1).FaceColor = [104 104 104]./255;
    b(2).FaceColor = [48 48 48]./255;
else
    b(1).FaceColor = [0 128 255]./255;
    b(2).FaceColor = [255 128 0]./255;
end
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
    lH = length(H);
    Htrim = [ H(1:min(lH,5)) zeros(1,5-lH) ];
    totT = [totT ; Htrim];  %#ok
    strleg{n} = strcat('Nnodes=',num2str(varList(n)));
end
b = bar(totT.');
xlabel('Number of Wi-Fi nodes','FontSize',14);
ylabel('Percentage (%)','FontSize',14);
title('Number of nodes per group - ABS1','FontSize',14);
lg = legend(strleg);
set(lg,'FontSize',14);
lc = findobj(lg, 'property', 'Color');
if conf.blackAndWhite
    b(4).FaceColor = [192 192 192]./255;
    b(3).FaceColor = [144 144 144]./255;
    b(2).FaceColor = [104 104 104]./255;
    b(1).FaceColor = [48 48 48]./255;
else
    b(4).FaceColor = [102 178 255]./255;
    b(3).FaceColor = [0 128 255]./255;
    b(2).FaceColor = [0 76 153]./255;
    b(1).FaceColor = [0 51 102]./255;
end
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
if conf.blackAndWhite
    for i = 2:11
    b(i).FaceColor = [26*(i-2) 26*(i-2) 26*(i-2)]./255;
    end
else
    % Brown scale
    b(2).FaceColor = [51 25 0]./255;
    b(3).FaceColor = [102 51 0]./255;
    b(4).FaceColor = [153 76 0]./255;
    b(5).FaceColor = [204 102 0]./255;
    b(6).FaceColor = [255 128 0]./255;
    b(7).FaceColor = [255 153 51]./255;
    b(8).FaceColor = [255 178 102]./255;
    b(9).FaceColor = [255 204 153]./255;
    b(10).FaceColor = [255 229 204]./255;
    b(11).FaceColor = [255 235 212]./255;
end
set(gca,'XTickLabel',[5 10 15 20] );
title('Latency comparison - LTE-U vs E-Fi','FontSize',15);
title('LTE-U vs E-Fi - Latency Analysis','FontSize',15);