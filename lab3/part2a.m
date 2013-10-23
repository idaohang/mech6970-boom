clear all; clc; close all;

part2c_load_data



x1=gNovatel0.zX_gNovatel0;
y1=gNovatel0.zY_gNovatel0;
z1=gNovatel0.zZ_gNovatel0;

x2=gNovatel1.zX_gNovatel1;
y2=gNovatel1.zY_gNovatel1;
z2=gNovatel1.zZ_gNovatel1;


diff=[x2-x1 y2-y1 z2-z1];

dist=sqrt(diff(:,1).*diff(:,1)+diff(:,2).*diff(:,2)+diff(:,3).*diff(:,3));

figure;
plot(time./1000,dist(1:616))
hold on;
plot(time./1000,1.905.*ones(1,length(time)),'r');
xlabel('Time')
ylabel('Position Difference')
legend('Measured Difference', 'True Difference')

nanmean(dist)
