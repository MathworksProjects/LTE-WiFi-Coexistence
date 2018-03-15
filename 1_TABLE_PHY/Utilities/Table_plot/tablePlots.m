clear; close all; clc;

% load('MCS0-ABS1.mat')
load('TABLE_PER_ABS0.mat')

% Old Variable Definition in TABLE
% x = Prx_WiFi_glob;
% y = SINR_WiFi_glob;
% z = 100-PER_WiFi_glob;

% New Variable Definition in TABLE
x = PRXLIST;
y = SINRLIST;
z = 100-PERLIST;

figure
scatter(x,y,10,z,'LineWidth',3.5);
% set(gca,'CLim',[0 1]);
xlabel('WiFi Received Power (dBm)','FontSize',14);
ylabel('WiFi SINR (dBm)','FontSize',14);
% zlabel('PER')
title('Packet Success Rate (%) of ABS1, MCS0','FontSize',14);

xlim([-95, -25])
ylim([-80, 60])
hold on

% for i = 1 : 2880
%     if(PER_WiFi_glob(i) <= 10)
%         plot(Prx_WiFi_glob(i), SINR_WiFi_glob(i), 'k*');
%     end
% end

% Plot limit lime for PSR threshold (tunnable)
t = -95:-25;
p = 9.8*ones(1, 71);
plot(t,p,'LineWidth',3);
annotation('textarrow',[0.2 0.2],[0.7 0.85],'String','$PSR \ge 90\%$', 'Interpreter','latex','FontSize',16,'FontWeight','bold');

colormap('gray')
colorbar
grid minor
box on

%% intrapolate
[xq, yq] = meshgrid(min(min(x), min(y)) : 1 : max(max(x), max(y)));
z1 = griddata(x,y,z,xq,yq,'cubic');
figure
plot3(x,y,z,'o')
hold on
surf(xq,yq,z1,'EdgeColor', 'none')
title('Cubic')
% legend('Sample Points','Interpolated Surface','Location','NorthWest')

xlim([-95, -25])
ylim([-80, 60])

grid minor
colormap('gray')
colorbar
xlabel('WiFi Received Power (dBm)')
ylabel('WiFi SINR (dBm)')
% zlabel('PER')
title('Packet Success Rate (%)')
