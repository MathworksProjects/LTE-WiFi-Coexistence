clear; close all; clc;

runnable = true;   %#ok % Always set to true if runnable

% Load the MAC Variables
run('../Configuration_MAC');

%% Throughput Experiment
PER_list        = (0:0.05:0.95);             % Do not Modify
conf.PSDULength = 600;                       % Default is 600
conf.Alpha      = 1;                         % Default is 1
conf.SaveTable  = true;                      % Default is False
conf.TableName  = 'TABLE_NEW_MAC_testing';   % Table Name
conf.Plotting   = true;                      % Controls Plotting
exp = 'PSDULength'; varList = (100:100:600);
% exp = 'Alpha'; varList = (0:0.1:1);
[~] = Throughput_experiment(exp, varList, PER_list, conf);

%% Closing Clause
runnable = false;

%% Function Throughput Experiment
function [out] = Throughput_experiment(experiment, varList, PER_list, conf)
PSDULength = conf.PSDULength;
alpha = conf.Alpha;
saveTable = conf.SaveTable;
plotting = conf.Plotting;

avTotalTime = zeros(length(PER_list),length(varList));
avTransmit = zeros(length(PER_list),length(varList));
avPktDropCnt = zeros(length(PER_list),length(varList));
for var = varList
    for PER = PER_list
        if(strcmp(experiment,'PSDULength'))
            PSDULength = var;
            [avTotalTime(PER==PER_list,var==varList), ...
                avTransmit(PER==PER_list,var==varList), ...
                avPktDropCnt(PER==PER_list,var==varList)] = ...
                throughputSimulation(PER, var, alpha);
            fprintf('%d %.2f - %.2f (ms) %.2f (ms) %.2f (%%)\n', ...
                var, PER, avTotalTime(PER==PER_list), ...
                avTransmit(PER==PER_list), avPktDropCnt(PER==PER_list));
        elseif(strcmp(experiment,'Alpha'))
            [avTotalTime(PER==PER_list,var==varList), ...
                avTransmit(PER==PER_list,var==varList), ...
                avPktDropCnt(PER==PER_list,var==varList)] = ...
                throughputSimulation(PER, PSDULength, var);
            fprintf('%.1f %.2f - %.2f (ms) %.2f (ms) %.2f (%%)\n', ...
                var, PER, avTotalTime(PER==PER_list), ...
                avTransmit(PER==PER_list), avPktDropCnt(PER==PER_list));
        end
    end
end

if plotting
    for var = varList
        figure(1); hold on; plot(PER_list,avTotalTime(:,var==varList));
        figure(2); hold on; plot(PER_list,PSDULength./avTotalTime(:,var==varList));
    end
    % Plot the Average time to Transmit in ms
    figure(1);
    plot(PER_list,1*ones(length(PER_list),1),'Color','k','LineStyle','--'); % 1ABS limit
    plot(PER_list,2*ones(length(PER_list),1),'Color','k','LineStyle','--'); % 2ABS limit
    plot(PER_list,3*ones(length(PER_list),1),'Color','k','LineStyle','--'); % 3ABS limit
    plot(PER_list,4*ones(length(PER_list),1),'Color','k','LineStyle','--'); % 4ABS limit
    plot(PER_list,5*ones(length(PER_list),1),'Color','k','LineStyle','--'); % 5ABS limit
    set(gca, 'XTick', PER_list);           % Change x-axis ticks
    set(gca, 'XTickLabel', 100.*PER_list); % Change x-axis ticks labels.
    hxlabel1 = xlabel('PER (%)');
    hylabel1 = ylabel('Average Time to transmit (ms)');
    if(strcmp(experiment,'PSDULength'))
        legendCell = cellstr(num2str(varList', 'PSDU(bits) = %-d'));
    elseif(strcmp(experiment,'Alpha'))
        legendCell = cellstr(num2str(varList', 'Alpha = %-.1f'));
    end
    hlegend1 = legend(legendCell,'Location','northwest');
    set(hxlabel1,'FontSize',14);
    set(hylabel1,'FontSize',14);
    set(hlegend1,'FontSize',14);
    grid minor;
    % Plot the Average Throughput in Kbps
    figure(2);
    set(gca, 'XTick', PER_list);           % Change x-axis ticks
    set(gca, 'XTickLabel', 100.*PER_list); % Change x-axis ticks labels.
    hxlabel1 = xlabel('PER (%)');
    hylabel1 = ylabel('Throughput (kbps)');
    if(strcmp(experiment,'PSDULength'))
        legendCell = cellstr(num2str(varList', 'PSDU(bits) = %-d'));
    elseif(strcmp(experiment,'Alpha'))
        legendCell = cellstr(num2str(varList', 'Alpha = %-.1f'));
    end
    hlegend1 = legend(legendCell,'Location','northwest');
    set(hxlabel1,'FontSize',14);
    set(hylabel1,'FontSize',14);
    set(hlegend1,'FontSize',14);
    grid minor;
    % Plot the Average Number of Transmissions in Pkts
    figure(3);
    plot(PER_list,avTransmit(:,1));
    set(gca, 'XTick', PER_list);           % Change x-axis ticks
    set(gca, 'XTickLabel', 100.*PER_list); % Change x-axis ticks labels.
    hxlabel2 = xlabel('PER (%)');
    hylabel2 = ylabel('Average Number of Transmissions (Pkts)');
    hlegend2 = legend('802.11ac Packet','Location','northwest');
    set(hxlabel2,'FontSize',14);
    set(hylabel2,'FontSize',14);
    set(hlegend2,'FontSize',14);
    grid minor;
    % Plot the Probability of Dropping a Packet in %
    figure(4);
    plot(PER_list,avPktDropCnt(:,1));
    set(gca, 'XTick', PER_list);           % Change x-axis ticks
    set(gca, 'XTickLabel', 100.*PER_list); % Change x-axis ticks labels.
    hxlabel3 = xlabel('PER (%)');
    hylabel3 = ylabel('Probability of Dropping a Packet');
    hlegend3 = legend('802.11ac Packet','Location','northwest');
    set(hxlabel3,'FontSize',14);
    set(hylabel3,'FontSize',14);
    set(hlegend3,'FontSize',14);
    grid minor;
end

% Create New Table for Throughput Calculations
if saveTable && strcmp(experiment,'PSDULength')
	lPERList  = length(PER_list);
    lPSDUList = length(varList);  % varList = PSDULength
    out.THPSDULIST = [];
    out.THTXTIMELIST = [];
    out.THTXNUMBLIST = [];
    out.THPDROPLIST = [];
    for idx = 1:lPSDUList
        out.THPSDULIST   = [out.THPSDULIST   repmat(varList(idx),1,lPERList)];
        out.THTXTIMELIST = [out.THTXTIMELIST avTotalTime(:,idx).'];
        out.THTXNUMBLIST = [out.THTXNUMBLIST avTransmit(:,idx).'];
        out.THPDROPLIST  = [out.THPDROPLIST  avPktDropCnt(:,idx).'];
    end
    out.THPERLIST = repmat(PER_list,1,lPSDUList);
elseif saveTable && ~strcmp(experiment,'PSDULength')
    disp('The Experiment should be set to PSDULength in order to create the table');
    return;
end
end

%% Function Throughput Simulation
function [avTotalTime, avTransmit, avPktDropCnt] = throughputSimulation(PER, PSDULength, alpha)

% 802.11n MAC parameters. Do not modify
DIFS = 34e-6;       % DIFS time
SIFS = 16e-6;       % SIFS time
slotTime = 9e-6;    % slot time
maxReTx = 7;        % Maximum retransmissions before dropping a Pkt

switch PSDULength
    case 100
        packetTime = 270e-6;
    case 200
        packetTime = 394e-6;
    case 300
        packetTime = 518e-6;
    case 400
        packetTime = 638e-6;
    case 500
        packetTime = 762e-6;
    case 600
        packetTime = 886e-6;
end

Niter = 100000;

rTransmit = zeros(Niter,1);
rPktDropCnt = zeros(Niter,1);
rTotalTime = zeros(Niter,1);

for iter = 1 : Niter
    transmitSuccess = false;
    contentionWindow = 15 * slotTime;
    packetDropCounter = 0;
    transmitCounter = 0;
    totalTime = 0;

    while ( ~transmitSuccess && transmitCounter <= maxReTx )
        backoffTime = alpha*(rand(1, 1) * contentionWindow);
        totalTime = totalTime + packetTime + alpha*DIFS + SIFS + backoffTime;
        randNumber = rand(1,1);
        if(randNumber >= PER)
            transmitSuccess = true;
        end
        contentionWindow = contentionWindow * 2;
        transmitCounter = transmitCounter + 1;
    end

    if(transmitCounter > maxReTx)
        % The packet is dropped. It doesn't have an effect on the time
        packetDropCounter = packetDropCounter + 1;
        rTotalTime(iter) = NaN;
        rTransmit(iter) = NaN;
    else
        rTotalTime(iter) = totalTime;
        rTransmit(iter) = transmitCounter;
    end
	rPktDropCnt(iter) = packetDropCounter;
end

rTotalTime1 = rTotalTime(~isnan(rTotalTime));    % Get ride of NaN
rTransmit1  = rTransmit(~isnan(rTransmit));      % Get ride of NaN

avTotalTime = mean(rTotalTime1)*1e3;             % In ms
avTransmit = mean(rTransmit1);                   % In Pkt
avPktDropCnt = 100*sum(rPktDropCnt)/Niter;       % In %
end