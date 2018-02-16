clear; close all; clc;

run('../Configuration_MAC');

tableName = strcat(conf.TablePath,conf.TableName);
load(tableName);

x = THPSDULIST;
y = THPERLIST;
z = THTXTIMELIST;

colorList = parula(3);

indx(1)    = length(ALPHAORG);            % Alpha = 1
color(1,:) = colorList(:, 1);
indx(2)    = ceil(length(ALPHAORG)/2);    % Alpha = 0.5
color(2,:) = colorList(:, 2);
indx(3)    = 1;                           % Alpha = 0
color(3,:) = colorList(:, 3);

figure; hold on;
stem3(x,y,z(indx(1),:),'s','fill','Color',color(1,:)); 
stem3(x,y,z(indx(2),:),'s','fill','Color',color(2,:));
stem3(x,y,z(indx(3),:),'s','fill','Color',color(3,:));
plot3(x,y,zeros(length(x),1), '.r');
xlabel('PSDU Length (bits)','FontSize',12);
ylabel('PER','FontSize',12);
zlabel('Average Time to Succesfully Tx (ms)','FontSize',12);
hleg = legend('Alpha=1','Alpha=0.5','Alpha=0');
set(hleg,'FontSize',12);
view([-31.9 20.4]);
grid minor;
title('Average Time to Succesfully transmit a Wi-Fi Packet');

%% Interpolate
[xq, yq] = meshgrid(unique(x),unique(y));
z1 = griddata(x,y,z(indx(1),:),xq,yq,'cubic');
z2 = griddata(x,y,z(indx(2),:),xq,yq,'cubic');
z3 = griddata(x,y,z(indx(3),:),xq,yq,'cubic');
figure; hold on;
plot3(x,y,z(indx(1),:),'s','MarkerSize',5,'MarkerFaceColor',color(1,:),'Color',color(1,:))
plot3(x,y,z(indx(2),:),'s','MarkerSize',5,'MarkerFaceColor',color(2,:),'Color',color(2,:))
plot3(x,y,z(indx(3),:),'s','MarkerSize',5,'MarkerFaceColor',color(3,:),'Color',color(3,:))
hSurf = surf(xq,yq,z1,'EdgeColor', 'none');
set(hSurf,'FaceColor',color(1,:),'FaceAlpha',0.4);
hSurf = surf(xq,yq,z2,'EdgeColor', 'none');
set(hSurf,'FaceColor',color(2,:),'FaceAlpha',0.6);
hSurf = surf(xq,yq,z3,'EdgeColor', 'none');
set(hSurf,'FaceColor',color(3,:),'FaceAlpha',0.8);

xlim([min(x), max(x)]);
ylim([min(y), max(y)]);

xlabel('PSDU Length (bits)','FontSize',12);
ylabel('PER','FontSize',12);
zlabel('Average Time to Succesfully Tx (ms)','FontSize',12);
hleg = legend('Alpha=1','Alpha=0.5','Alpha=0');
set(hleg,'FontSize',12);
view([-31.9 20.4]);
grid minor;
title('Average Time to Succesfully transmit a Wi-Fi Packet');