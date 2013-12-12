%% Numerically Integrate Crossbow

clear all; close all; clc
load workspace_after_time_sync

[b,a]=butter(4, .005, 'low');


for k = 1:4
  XBow_accY_filt(:,k) = filtfilt(b,a,XBow_accY(:,k));
end

% remove gravity
XBow_accY_filt = XBow_accY_filt - mean(XBow_accY_filt(1:1000,1));

plot(XBow_accY(:,1)); hold on
plot(XBow_accY_filt(:,1), 'g')

dTOW = mean(mean(diff(TOW(:,1:4))));
for ch = 1:4
  XBow_velY(:,ch) = cumtrapz(XBow_accY_filt(:,ch)) * dTOW;
end
XBow_velY = XBow_velY(1:end-1,:);
XBow_velY = [zeros(1,4); XBow_velY];

figure;
plot(XBow_velY(:,1))


%% calc doppler

userpos_lla = [gNovatel.zLat(1), gNovatel.zLong(1), gNovatel.zHeight(1)];

stop_idx =  50886; % data is bad after this

% interpolate GPS course
raw_course = unwrap(gNovatel.zCourse/360*2*pi);
imu_dopp = zeros(length(TOW),4);
% course = zeros(length(TOW),4);

for k = 1
  course = interp1(gNovatel.zGPSSeconds/1000, raw_course, TOW(:,k));
  course = rem(course(:,k),2*pi); 
  for kk = 1:length(TOW)
    imu_dopp(kk,k) = calcDoppler(userpos_lla,course(kk), XBow_velY(kk,k), svvel_nordnav{kk,k}, svpos_nordnav{kk,k});
  end
  
end

f = figure;
plot(TOW(:,1),imu_dopp(:,1))
grid on
title('Doppler resultant from integrated IMU');
legend('32','16','20','31')
ylabel('f_{Dopp} (Hz)')
xlabel('GPS Time (sec)');
saveas(f, 'Dopplers_over_time.png')
