function [ medianCQI, medianCQI_ABS, throughput_glob] = lte_report(totalTBS, totalCRC, CQIReport, CQIReport_ABS, numFrames)
    %% THROUGPUT LTE
    tputUEReported = sum(totalTBS.*(1-totalCRC));

    medianCQI = ceil(median(CQIReport));
    medianCQI_ABS = ceil(median(CQIReport_ABS));

    throughput_glob = (tputUEReported/(numFrames*1e-2))*1e-6;
end