clear all; close all; clc

% Settings 

% alog_mat = 'run_south_alog.mat';
% output_mat = 'run_south_alog_interp.mat';
alog_mat = 'run_north_alog.mat';
output_mat = 'run_north_alog_interp.mat';

load(alog_mat)

% make time vectors agree
start_time = min(gXbow440.time(find(gXbow440.time>min(gNovatel.time))));
stop_time = max(gXbow440.time(find(gXbow440.time<max(gNovatel.time))));
start_time_valid_idx = find(gXbow440.time>=start_time);
stop_time_valid_idx = find(gXbow440.time<=stop_time);
valid_idx = intersect(start_time_valid_idx,stop_time_valid_idx);
gXbow440.time = gXbow440.time(valid_idx);
gXbow440.zGyroX = gXbow440.zGyroX(valid_idx);
gXbow440.zGyroY = gXbow440.zGyroY(valid_idx);
gXbow440.zGyroZ = gXbow440.zGyroZ(valid_idx);
gXbow440.zAccelX = gXbow440.zAccelX(valid_idx);
gXbow440.zAccelY = gXbow440.zAccelY(valid_idx);
gXbow440.zAccelZ = gXbow440.zAccelZ(valid_idx);


% clean up nans
valid_idx = find(~isnan(gNovatel.zVertVel));
gNovatel.time = gNovatel.time(valid_idx);
gNovatel.zLat = gNovatel.zLat(valid_idx);
gNovatel.zLong = gNovatel.zLong(valid_idx);
gNovatel.zHeight = gNovatel.zHeight(valid_idx);
gNovatel.zHorizSpeed = gNovatel.zHorizSpeed(valid_idx);
gNovatel.zVertVel = gNovatel.zVertVel(valid_idx);
gNovatel.zCourse = gNovatel.zCourse(valid_idx);
gNovatel.zGPSWeek = gNovatel.zGPSWeek(valid_idx);
gNovatel.zGPSSeconds = gNovatel.zGPSSeconds(valid_idx);


% interpolation
gXbow440.gps_time = interp1( gNovatel.time, gNovatel.zGPSSeconds, gXbow440.time , 'cubic');

%

save(output_mat, 'gNovatel','gXbow440')