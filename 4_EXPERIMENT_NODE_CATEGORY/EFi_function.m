function [avThpt,avPSR,nodes,PSR,PSR_new,posNodes,outFeasExp]= EFi_function(conf)
% Copy Variables to Local Variables
Radius      = conf.Radius;
Cover_angle = conf.Cover_angle;
posAP       = conf.posAP;
posBS       = conf.posBS;
Ptx_WiFi    = conf.Ptx_WiFi;
Ptx_LTE     = conf.Ptx_LTE;
Pnoise      = conf.Pnoise;
Nnodes      = conf.Nnodes;
PSR_th      = conf.PSR_th;
K           = conf.K;
PSDULength  = conf.PSDULength;
ABSCfg      = conf.ABSCfg;
Alpha       = conf.Alpha;
Criteria    = conf.Criteria;
NS3Bridge   = conf.NS3Bridge;
NS3FileName = conf.NS3FileName;
PSR_minth   = conf.PSR_minth;
PSR_maxth   = conf.PSR_maxth;

outFeasExp.cand = [];
outFeasExp.final = [];

% Load the PER Table if not done before
if(~isfield(conf,'PERLIST'))
    if(ABSCfg == 0);     load('TABLE_PER_ABS0'); 
    elseif(ABSCfg == 1); load('TABLE_PER_ABS1');
    end
    PERLIST  = PERLIST./100;
else
    PERLIST  = conf.PERLIST;
    PRXLIST  = conf.PRXLIST;
    SINRLIST = conf.SINRLIST;
end

% Load the MAC Table if not done before
if(~isfield(conf,'THPERLIST'))
%     load('TABLE_NEW_MAC');
    load('TABLE_NEW_MAC_testing');
	% Times for Alpha=variable (Modified) (AP->GO and GO->P2P)
    [~,idxAPtoGO]  = min(abs(repmat(Alpha,size(ALPHAORG)) - ALPHAORG));
    [~,idxGOtoP2P] = min(abs(repmat(1-Alpha,size(ALPHAORG)) - ALPHAORG));
    THTXTIMELIST_APtoGO = THTXTIMELIST(idxAPtoGO,:);   %#ok % AP->GO
    THTXTIMELIST_GOtoP2P = THTXTIMELIST(idxGOtoP2P,:);      % GO->P2P
	% Times for Alpha=1 (Not Modified) (GO, CSZ and CNSZ)
    THTXTIMELIST = THTXTIMELIST(end,:);                     % AP -> Rest
else
    THPERLIST = conf.THPERLIST;
    THPSDULIST = conf.THPSDULIST;
    THTXTIMELIST = conf.THTXTIMELIST.*1e-3;
    THTXTIMELIST_APtoGO = conf.THTXTIMELIST_APtoGO.*1e-3;
    THTXTIMELIST_GOtoP2P = conf.THTXTIMELIST_GOtoP2P.*1e-3;
end

% For Throughput calculation purposes, in case the TABLE_NEW_MAC is not
% used. This serves as an approximation, and it is a lower bound of the
% more realistic calculation of the Throughput.
ABSDuration = 1e-3;                     % 1ms
if(PSDULength<200);                              NumberOfPacketPerABS = 3;
elseif( (PSDULength>=200) && (PSDULength<300) ); NumberOfPacketPerABS = 2;
elseif( (PSDULength>=300) && (PSDULength<700) ); NumberOfPacketPerABS = 1;
else;                                            NumberOfPacketPerABS = 0;
end
rate = ABSDuration/NumberOfPacketPerABS;

%% Random locations of Nnodes within the Coverage Area define by the Radius
minR = 0; maxR = Radius;
minAg = -Cover_angle;
maxAg = Cover_angle;
radius = (maxR-minR).*rand(Nnodes,1) + minR;
angle = (maxAg-minAg).*rand(Nnodes,1) + minAg;
posNodes = [radius.*cos(angle) radius.*sin(angle)];

%% Calculate Respective distance to the AP and to the BS
D_WiFi = sqrt((posNodes(:,1)-posAP(1)*ones(Nnodes,1)).^2 ...
    + (posNodes(:,2)-posAP(2)*ones(Nnodes,1)).^2);
D_LTE = sqrt((posNodes(:,1)-posBS(1)*ones(Nnodes,1)).^2 ...
    + (posNodes(:,2)-posBS(2)*ones(Nnodes,1)).^2);

%% Calculate the WiFi and LTE RX Power and the SINR
Pf = 0;     % Floor Penetration loss factor for a regular office in the 5GHz band
            % When transmitting and receiving in the same floor.
N = 31;     % Distance Power Loss Coeficient for an office Area in the 5GHz band
loss_WiFi = -(20*log10(5.2e3) - 28 + N*log10(D_WiFi) + Pf);
PRX_WiFi = Ptx_WiFi + loss_WiFi;
PRX_WiFi_lin = 10.^((PRX_WiFi-30)./10);

loss_LTE = -(20*log10(5.2e3) - 28 + N*log10(D_LTE) + Pf);
PRX_LTE = Ptx_LTE + loss_LTE;
PRX_LTE_lin = 10.^((PRX_LTE-30)./10);

Pnoise_lin = 10.^((Pnoise-30)./10);
SINR_lin = PRX_WiFi_lin./(PRX_LTE_lin+Pnoise_lin);
SINR = 10*log10(SINR_lin);

%% Map the Results to Obtain the PSR as a function of the RX Power and SINR
PER = [];
for n = 1:Nnodes
    dist = sqrt( (ones(1,size(PRXLIST,2))*PRX_WiFi(n) - PRXLIST).^2 + ...
           (ones(1,size(SINRLIST,2))*SINR(n) - SINRLIST).^2);
    [~,ltn] = min(dist);
    PER = [PER PERLIST(ltn)];                       %#ok<AGROW>
end

% Convert PER from table into PSR for better visualization
PSR = 1 - PER;  

% Criteria on how to elect the Group Owner and Wi-Fi Direct Candidates.
% This assignation is not final, since the Hungarian algorithm will
% determine wheather or not to change their role.
GOCList = find(PSR>=PSR_th);
P2PCList = find(PSR<PSR_th);

if((~isempty(GOCList) && ~isempty(P2PCList)))
    %% PSR through Candidate Relays for all P2P
    posP2PCNodes = [posNodes(P2PCList,1) posNodes(P2PCList,2)];
    posGOCCNodes = [posNodes(GOCList,1) posNodes(GOCList,2)];
    PRX_LTE1 = PRX_LTE(P2PCList);
    PRX_LTE_lin1 = 10.^((PRX_LTE1-30)./10);
    PSR_new = [];
    PSR_GOtoP2P = [];
    for p2pc = 1:length(P2PCList)
        dist_W = sqrt((posGOCCNodes(:,1)-posP2PCNodes(p2pc,1)*ones(length(GOCList),1)).^2 + ...
                (posGOCCNodes(:,2)-posP2PCNodes(p2pc,2)*ones(length(GOCList),1)).^2);
        loss_WiFi1 = -(20*log10(5.2e3) - 28 + N*log10(dist_W) + Pf);
        PRX_WiFi1 = Ptx_WiFi + loss_WiFi1;
        PRX_WiFi_lin1 = 10.^((PRX_WiFi1-30)./10);

        SINR_lin1 = PRX_WiFi_lin1./(PRX_LTE_lin1(p2pc)+Pnoise_lin);
        SINR1 = 10*log10(SINR_lin1);

        PER1 = [];
        for n = 1:length(GOCList)
            dist = sqrt( (ones(1,size(PRXLIST,2))*PRX_WiFi1(n) - PRXLIST).^2 + ...
                (ones(1,size(SINRLIST,2))*SINR1(n) - SINRLIST).^2);
            [~,ltn] = min(dist);
            PER1 = [PER1 PERLIST(ltn)];            %#ok<AGROW>
        end
        PSR1 = (1-PER1);

        % Evaluate Equation in the Paper to define the candidates
        if(strcmp(Criteria,'PER_or'))        % Criteria 1
            PSR1 = PSR1.*PSR(GOCList);
            noncandidates = find(PSR1<=PSR(P2PCList(p2pc)));
        elseif(strcmp(Criteria,'Th_or'))     % Criteria 2
            PSR2 = 1./(1./PSR1 + 1./PSR(GOCList));
            noncandidates = find(PSR2<=PSR(P2PCList(p2pc)));
            PSR1 = PSR1.*PSR(GOCList);
        end

        PSR1(noncandidates) = 0;
        PSR_new = [PSR_new ; PSR1];                %#ok<AGROW>
        PSR_GOtoP2P = [PSR_GOtoP2P ; (1-PER1)];    %#ok<AGROW>
    end
    
    list = find( (PSR>=PSR_minth) & (PSR<=PSR_maxth) );
    outFeasExp.cand = [];
    for l = list
        x = find(l==P2PCList);
        list2 = find(PSR_new(x,:)~=0);
        if ~isempty(list2)
            for k = list2
                p = GOCList(k);

                thptP2P     = calculateThroughputClient(conf, PSR(l), THPERLIST, THPSDULIST, THTXTIMELIST, rate);
                thptP2P_new = calculateThroughputWiFiDirect(conf, PSR(p), PSR_GOtoP2P(x,k), THPERLIST, THPSDULIST, THTXTIMELIST_APtoGO, THTXTIMELIST_GOtoP2P, rate);

                outFeasExp.cand = [outFeasExp.cand ; PSR(l) PSR(p) PSR_GOtoP2P(x,k) thptP2P*1e-3 thptP2P_new*1e-3]; %#ok<AGROW>
            end
        end
    end

    %% Cost Matrix as the Input for the Hungarian Algorithm
    W = [PSR_new PSR(P2PCList).'];
    W_GOtoP2P = [PSR_GOtoP2P PSR(P2PCList).'];

    %% Initialization of the Modified Hungarian Algorithm
    M = length(GOCList); N = length(P2PCList);
    P = ceil(ones(size(W,1),size(W,2)) - (ones(size(W,1),size(W,2)) - W).^50);
    t = ones(1,size(P,1))*P;
    dummy_relays = t;

    % Dummy Relays insertion
    W_Mod = []; pVec = [];
    for col = 1:size(W,2)-1    %The last one belongs to the AP
        W1 = repmat(W(:,col),1, min(dummy_relays(col),K));
        W_Mod = [W_Mod W1];                                   %#ok<AGROW>
        pV1 = repmat(GOCList(col), 1, min(dummy_relays(col),K));
        pVec = [pVec pV1];                                    %#ok<AGROW>
    end

    % Dummy AP insertion
    W2 = repmat(W(:,size(W,2)), 1, Nnodes-M);
    W_Mod = [W_Mod W2];
    pV1 = repmat(-1, 1, Nnodes-M);
    pVec = [pVec pV1];

    %% Hungarian Algorithm
    W_Mod_reverse = ones(size(W_Mod,1),size(W_Mod,2)) - W_Mod;
    [ Assign , ~ ] = munkres(W_Mod_reverse);
    P2PGOList = pVec(Assign);           % (X,f(X)) = P2P X assigned to GOC f(X)

    %% Final PSR in the Network
    PSR_new = PSR;
    PSR_GOtoP2P = PSR;

    % PSR for the P2PC
    PSR_new_p2p = [];
    PSR_new_p2p1 = [];
    for p2pc = 1:length(P2PCList)
        if(P2PGOList(p2pc)~=-1)
            loc_M = find(GOCList==P2PGOList(p2pc));
            PSR_new_p2p = [PSR_new_p2p W(p2pc,loc_M)];           %#ok<AGROW>
            PSR_new_p2p1 = [PSR_new_p2p1 W_GOtoP2P(p2pc,loc_M)]; %#ok<AGROW>
        else
            PSR_new_p2p = [PSR_new_p2p W(p2pc,end)];             %#ok<AGROW>
            PSR_new_p2p1 = [PSR_new_p2p1 W_GOtoP2P(p2pc,end)];   %#ok<AGROW>
        end
    end
    PSR_new(P2PCList) = PSR_new_p2p;
    PSR_GOtoP2P(P2PCList) = PSR_new_p2p1;

    %% Final Node Categorization
    % P2PC Relaying
	relay = (1:Nnodes);
    relay(P2PCList) = P2PGOList;

    % Final Node Grouping
    relay_final = relay;
    loc = find(relay==-1);
    relay_final(loc) = loc;

    % Clients N.S.Z
    loc = find(relay(P2PCList)==-1);
    Clients_NSZ = P2PCList(loc);
    % Clients S.Z VS GO
    Clients_SZ = []; GO = []; Psi = [];
    for goc = GOCList
        Psi_List = find(relay==goc);
        Psi_List(Psi_List==goc) = [];
        if(isempty(Psi_List))
            % Clients S.Z
            Clients_SZ = [Clients_SZ goc];    %#ok<AGROW>
        else
            % GO
            GO  = [GO goc];                   %#ok<AGROW>
            Psi = [Psi length(Psi_List)];     %#ok<AGROW>
        end
    end
    % P2P Nodes
    P2P = setdiff((1:Nnodes), [Clients_SZ Clients_NSZ GO]);
else
    %% Final Node Categorization
    PSR_new = PSR;

    GO = []; P2P = []; Psi = [];
    Clients_SZ = GOCList;
    Clients_NSZ = P2PCList;
    P2PGOList = P2PCList;

    relay = (1:Nnodes);
    relay(P2PCList) = P2PGOList;
    relay_final = relay;
end

% Node Categorization
nodes.GO = GO;
nodes.P2P = P2P;
nodes.Clients_SZ = Clients_SZ;
nodes.Clients_NSZ = Clients_NSZ;

% Average PSR
avPSR.old = mean(PSR);
avPSR.new = mean(PSR_new);
% avPSR.old = mean(PSR(P2P));
% avPSR.new = mean(PSR_new(P2P));

%% Throughput Calculation

thptGO = zeros(length(GO),1);
thptP2P = zeros(length(P2P),1);
thptP2P_new = zeros(length(P2P),1);
thptCSZ = zeros(length(Clients_SZ),1);
thptCNSZ = zeros(length(Clients_NSZ),1);
for go = GO
    thptGO(go==GO) = calculateThroughputClient(conf, ...
        PSR(go), THPERLIST, THPSDULIST, THTXTIMELIST, rate);
end
for p2p = P2P    
    thptP2P(p2p==P2P) = calculateThroughputClient(conf, ...
        PSR(p2p), THPERLIST, THPSDULIST, THTXTIMELIST, rate);
    thptP2P_new(p2p==P2P) = calculateThroughputWiFiDirect(conf, ...
        PSR(relay(p2p)), PSR_GOtoP2P(p2p), THPERLIST, THPSDULIST, ...
        THTXTIMELIST_APtoGO, THTXTIMELIST_GOtoP2P, rate);
end
for csz = Clients_SZ
    thptCSZ(csz==Clients_SZ) = calculateThroughputClient(conf, ...
        PSR(csz), THPERLIST, THPSDULIST, THTXTIMELIST, rate);
end
for cnsz = Clients_NSZ
    thptCNSZ(cnsz==Clients_NSZ) = calculateThroughputClient(conf, ...
        PSR(cnsz), THPERLIST, THPSDULIST, THTXTIMELIST, rate);
end

avThpt.GO = mean(thptGO)*1e-3;
avThpt.CSZ = mean(thptCSZ)*1e-3;
avThpt.CNSZ = mean(thptCNSZ)*1e-3;
avThpt.P2P = mean(thptP2P)*1e-3;
avThpt.P2P_new = mean(thptP2P_new)*1e-3;

list = find( (PSR>=PSR_minth) & (PSR<=PSR_maxth) );
outFeasExp.final = [];
if ~isempty(P2P)
    for p2p = list
        r = relay_final(p2p);

        thptP2P     = calculateThroughputClient(conf, PSR(p2p), THPERLIST, THPSDULIST, THTXTIMELIST, rate);
        thptP2P_new = calculateThroughputWiFiDirect(conf, PSR(r), PSR_GOtoP2P(p2p), THPERLIST, THPSDULIST, THTXTIMELIST_APtoGO, THTXTIMELIST_GOtoP2P, rate);

        outFeasExp.final = [outFeasExp.final ; PSR(p2p) PSR(r) PSR_GOtoP2P(p2p) thptP2P*1e-3 thptP2P_new*1e-3];
    end
end

%% Resource Assignation per Wi-Fi Group
loadFactor = zeros(length(GO)+1,1)';
for golist = 1:length(GO)
	list = find(relay_final==GO(golist));
    loadFactor(golist) = sum(PSR_new(list));
end

list = setdiff((1:Nnodes), [GO]);
loadFactor(end) = sum(PSR_new(list));
loadFactor = 100*loadFactor./sum(loadFactor);

%% Generate the Bridge File for NS-3
if NS3Bridge
    FID = fopen(NS3FileName, 'a');
    nodeType = zeros(Nnodes,1);
    nodeType(GO) = 1;
    nodeType(P2P) = 2;
    nodeType(Clients_SZ) = 3;
    nodeType(Clients_NSZ) = 3;
    for nodeID = 1:Nnodes
        fprintf(FID,'%2.0f %2.0f %4.1f %5.1f %5.1f %5.1f %5.1f %2.0f ',...
            nodeID,nodeType(nodeID),...
            posNodes(nodeID,1),posNodes(nodeID,2),...
            100*PSR(nodeID),100*PSR_GOtoP2P(nodeID),100*PSR_new(nodeID),...
            relay_final(nodeID));
    end
    fprintf(FID,'\n');

    for goIndx = 1:length(GO)+1
        fprintf(FID,'%2.2f ',loadFactor(goIndx));
    end
    fprintf(FID,'\n');
end

end

% Throughput calculation for Group Owners (GO), Clients in the Safe Zone
% (CSZ) and Clients not in the Safe Zone (CNSZ). The devices communicate
% directly with the AP, with no relying mechanisms.
function throughput = calculateThroughputClient(conf, PSR, THPERLIST, ...
    THPSDULIST, THTXTIMELIST, rate)

    PSDULength_mod = conf.PSDULength*(2^conf.MCS);

    % Old Method (Approximation, without taking into accound MAC Procedure)
%     time = rate * 1/PSR;
%     throughput = PSDULength_mod / time;
    
    % New Method using TABLE_NEW_MAC
    dist = sqrt( (ones(1,size(THPERLIST,2))*(1-PSR) - THPERLIST).^2 + ...
            (ones(1,size(THPERLIST,2))*PSDULength_mod - THPSDULIST).^2);
    [~,ltn] = min(dist);
    throughput = PSDULength_mod / THTXTIMELIST(ltn);
end

% Throughput calculation for Wi-Fi Direct Clients (P2P). The devices perform a
% 1-hop communication through a relay (Group Owner or GO). The throughput
% is determined based on the time to transmit from the AP to the GO plus
% the time to transmit from the GO to the P2P. In turn, these times are
% impacted by the alpha factor, which controls the modified Backoff
% mechanism.
function throughput = calculateThroughputWiFiDirect(conf, PSRAPtoGO, PSRGOtoP2P, ...
    THPERLIST, THPSDULIST, THTXTIMELIST_APtoGO, THTXTIMELIST_GOtoP2P, rate)
    
    PSDULength_mod = conf.PSDULength*(2^conf.MCS);
    
    % Old Method (Approximation, without taking into accound MAC Procedure)
%     time = rate * ( 1/PSRAPtoGO + 1/PSRGOtoP2P );
%     PSRGOtoP2P = PSDULength_mod / time;
    
    % New Method using TABLE_NEW_MAC
    dist = sqrt( (ones(1,size(THPERLIST,2))*(1-PSRAPtoGO) - THPERLIST).^2 + ...
                (ones(1,size(THPERLIST,2))*PSDULength_mod - THPSDULIST).^2);
    [~,ltn1] = min(dist);
    dist = sqrt( (ones(1,size(THPERLIST,2))*(1-PSRGOtoP2P) - THPERLIST).^2 + ...
                (ones(1,size(THPERLIST,2))*PSDULength_mod - THPSDULIST).^2);
    [~,ltn2] = min(dist);
    throughput = PSDULength_mod / (THTXTIMELIST_APtoGO(ltn1) + THTXTIMELIST_GOtoP2P(ltn2));
end