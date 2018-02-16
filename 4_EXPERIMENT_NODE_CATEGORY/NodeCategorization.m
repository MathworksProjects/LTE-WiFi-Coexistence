clear; close all; clc;

%% ------------------------------ VARIABLES -----------------------------%%
runnable = true;  %#ok % Always set to true if runnable

% Load the MAC Variables
run('../Configuration_MAC');

%% Node Categorization Experiment
conf.Niter = 200;        % Keep it above 200 iterations for meaningful results
PSRthList = [0.9 0.6];   % PSR threshold list to evaluate in the simulations
Categorization_experiment(conf,PSRthList);

%% Closing Clause
runnable = false;

%% Function Node Categorization
function [] = Categorization_experiment(conf,PSRthList)
    Niter = conf.Niter;
    alpha = conf.Alpha;
    
	% Load PSR/PER Tables (PRX,SINR) -> PSR
    if(conf.ABSCfg == 0);      tablePSRName = strcat(conf.TablePath,'TABLE_PER_ABS0');
    elseif(conf.ABSCfg == 1)
        if(conf.MCS == 0);     tablePSRName = strcat(conf.TablePath,'TABLE_PER_ABS1');
        elseif(conf.MCS == 1); tablePSRName = strcat(conf.TablePath,'TABLE_MCS1_ABS1');
        elseif(conf.MCS == 2); tablePSRName = strcat(conf.TablePath,'TABLE_MCS2_ABS1');
        elseif(conf.MCS == 3); tablePSRName = strcat(conf.TablePath,'TABLE_MCS3_ABS1');
        end
    end
    load(tablePSRName,conf.varPSRList{:});
    conf.PRXLIST  = PRXLIST;
    conf.SINRLIST = SINRLIST;
    conf.PERLIST  = PERLIST./100;

    % Load MAC Tables (PSR) -> Average time to Transmit
    tableMACName = strcat(conf.TablePath,conf.TableName);
    load(tableMACName,conf.varMACList{:});
    conf.THPERLIST = THPERLIST;
    conf.THPSDULIST = THPSDULIST;
    [~,idxAPtoGO]  = min(abs(repmat(alpha,size(ALPHAORG)) - ALPHAORG));
    [~,idxGOtoP2P] = min(abs(repmat(1-alpha,size(ALPHAORG)) - ALPHAORG));
    conf.THTXTIMELIST_APtoGO = THTXTIMELIST(idxAPtoGO,:);   %#ok % AP->GO
    conf.THTXTIMELIST_GOtoP2P = THTXTIMELIST(idxGOtoP2P,:);      % GO->P2P
    conf.THTXTIMELIST = THTXTIMELIST(end,:);                     % AP->P2P (alpha=1)
    
    for idx = 1:length(PSRthList)
	conf.PSR_th = PSRthList(idx);
    conf.PER_th = 1-conf.PSR_th;
    
    colorGO       = [0 0 0]./255;
    colorCSZ      = [144 144 144]./255;
    colorCNSZ     = [192 192 192]./255;
    colorP2P_old1 = [72 72 72]./255;
    colorP2P_old2 = [255,165,0]./255;
    colorP2P_new  = [104 104 104]./255;
    for iter = 1:Niter
        [~,~,nodes,PSR,PSR_new,posNodes] = EFi_function(conf);

        % Plot Results Characterization - PSR OLD
        figure(1); subplot(ceil(length(PSRthList)/2),2,idx); hold on;
        stem3(posNodes(nodes.GO,1),posNodes(nodes.GO,2),PSR(nodes.GO),'s','LineWidth',2,'Color',colorGO)
        stem3(posNodes(nodes.Clients_SZ,1),posNodes(nodes.Clients_SZ,2),PSR(nodes.Clients_SZ),'^','fill','Color',colorCSZ)
        stem3(posNodes(nodes.P2P,1),posNodes(nodes.P2P,2),PSR(nodes.P2P),'o','LineWidth',2,'Color',colorP2P_old1)
        stem3(posNodes(nodes.Clients_NSZ,1),posNodes(nodes.Clients_NSZ,2),PSR(nodes.Clients_NSZ),'*','Color',colorCNSZ)
        % Plot Results Characterization - PSR NEW
        figure(2); subplot(ceil(length(PSRthList)/2),2,idx); hold on;
        stem3(posNodes(nodes.GO,1),posNodes(nodes.GO,2),100.*PSR_new(nodes.GO),'s','fill','Color',colorGO)
        stem3(posNodes(nodes.Clients_SZ,1),posNodes(nodes.Clients_SZ,2),100.*PSR_new(nodes.Clients_SZ),'s','fill','Color',colorCSZ)
        stem3(posNodes(nodes.P2P,1),posNodes(nodes.P2P,2),100.*PSR_new(nodes.P2P),'s','fill','Color',colorP2P_new)
        stem3(posNodes(nodes.P2P,1),posNodes(nodes.P2P,2),100.*PSR(nodes.P2P),'s','fill','Color',colorP2P_old2)
        stem3(posNodes(nodes.Clients_NSZ,1),posNodes(nodes.Clients_NSZ,2),100.*PSR_new(nodes.Clients_NSZ),'s','fill','Color',colorCNSZ)
    end
    
    % Plot Coverage area (Initial Coverage Area)
    if(conf.Cover_angle==pi/2)
        axis([-2,conf.Radius,-conf.Radius,conf.Radius])
        x=(0:.01:conf.Radius);
    else
        axis([-conf.Radius,conf.Radius,-conf.Radius,conf.Radius])  
        x=(-conf.Radius:.01:conf.Radius);
    end
    y=[sqrt(conf.Radius.^2-x.^2) ; -sqrt(conf.Radius.^2-x.^2)];
    figure(1); subplot(ceil(length(PSRthList)/2),2,idx); hold on;
    plot(x,y,'Color','k'); 
    hleg = legend('Group Owner','Client S.Z.','P2P N.S.Z.','Client N.S.Z.');
    set(hleg,'Location','NorthEast');
    grid minor;
    figure(2); subplot(ceil(length(PSRthList)/2),2,idx); hold on;
    plot(x,y,'Color','k');
    grid minor;

    % Plot Safe Zone (modified Coverage Area) (2D or 3D)
    x = (0:2:conf.Radius);
    y = (-conf.Radius:2:conf.Radius);
    [X,Y] = meshgrid(x,y);                % Generate X and Y data
    Z = 100.*conf.PSR_th.*ones(size(X));  % Generate Z data
    figure(1); subplot(ceil(length(PSRthList)/2),2,idx); hold on;
    fig2H = copy(gcf);
    s = surf(X, Y, Z); set(s,'FaceAlpha',0.7); % Plot the surface
    view([17.76,26.00]);
    % view([0,0]);
    title('PSR (%)for Wi-Fi')
    hleg = legend('Group Owner','Client S.Z.','P2P N.S.Z.','Client N.S.Z.');
    set(hleg,'Location','NorthEast');
    figure(2); subplot(ceil(length(PSRthList)/2),2,idx); hold on;
    s = surf(X, Y, Z); % Plot the surface
    set(s,'FaceAlpha',0.7);
    view([17.76,26.00]);
    % view([0,0]);
    title('PSR (%) for E-Fi')
    hleg = legend('Group Owner','Client S.Z.','P2P N.S.Z. NEW','P2P N.S.Z. OLD','Client N.S.Z.');
    set(hleg,'Location','NorthEast');
    end
end