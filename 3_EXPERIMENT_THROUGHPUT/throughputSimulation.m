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