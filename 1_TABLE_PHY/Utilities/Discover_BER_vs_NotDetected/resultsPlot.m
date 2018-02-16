clc;clear;close all;
load('results_MCS0123_ALL.mat');

list = [2 4 8 14 25 40 54];
SINRlist = [50 40 30 20 10 0 -10];
MCS0 = [SINR_WiFi_glob(1:79); not_detected_glob(1:79); BER_glob(1:79)]';
MCS1 = [SINR_WiFi_glob(80:158); not_detected_glob(80:158); BER_glob(80:158)]';

MCS = [MCS0 MCS1];

figure
yyaxis left
hold on
plot(SINRlist,MCS0(list, 2),'-o');
plot(SINRlist,MCS1(list, 2),'-o');
ylabel('Not detect (%)')

yyaxis right
plot(SINRlist,MCS0(list, 3),'-o');
plot(SINRlist,MCS1(list, 3),'-o');
ylabel('Bit Error Rate (%)')

legend('MCS0, ND', 'MCS0, ND', 'MCS1, BER', 'MCS1, BER');

title('Not detected vs Bit Error Rate for ABS-1')
xlabel('SINR, dB');

grid minor;