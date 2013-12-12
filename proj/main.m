%% main.m Primary File for Final Project
% 
% 
clear all; close all; clc
tic

%% Settings / Constants

% run = 'south';
run = 'north';

% filter_order = 1;
% filter_order = 2;
% filter_order = 3;

data_dir = ['..' filesep 'data' filesep 'final_proj_data' filesep];

stop_time = 119; % how long to calculate in sec

c = 299792458; % speed of light (m/s)


%% load data

load(['run_' run '_alog_interp.mat']);
nordnav_fid = fopen([data_dir 'run_' run '.sim'], 'r');


%% Get SV Geometry Data
% Calculate pos, vel for all SVs for which we have ephemeris on the novatel

% % Get Ephemeris data
% use the Novatel and find it for each GPS time epoch
svs_have_ephem = [1 6 14 16 20 22 23 25 29 31 32];
nsv_have_ephem = length(svs_have_ephem);

% will store ephemeris as it came from the Novatel
ephem_novatel = zeros(30,nsv_have_ephem);
% % will store ephemeris in the format needed for the function `calc_sv_pos`
% ephem_gavlab = zeros(21,nsv_have_ephem);

% get data from Novatel
data_count = 1;
for k = 1:length(gNovatel.zPsrL1)
  % didn't log novatel raw data as fast as I thought.. oops
  if isempty(gNovatel.zPsrL1{k}), continue; end
  % calculate the transit time using Novatel's PSRs
  for ch = 1:32
    transit_time(data_count,ch) = gNovatel.zPsrL1{k}(ch+1) / c;
  end
  % get range of GPS times for which to calculate the SV Positions
  svpos_gpstimes(data_count) = gNovatel.zPsrL1{k}(34)/100;
  data_count = data_count+1;
end
clear data_count k
  
% sv positions at those GPS times
%   - rows are each time, cols are each sv
svpos = cell(length(svpos_gpstimes), nsv_have_ephem); 
svvel = svpos;
% SV clock corrections, dimension correspond to `svpos`
svclkcorr = zeros(size(svpos));


for ch = 1:nsv_have_ephem
  
  eph_ = getfield(gNovatel, ['zEphem' num2str(svs_have_ephem(ch))]);
  ephem_novatel(1:length(eph_{end}),ch) = eph_{end};
  %   ephem_gavlab(:,ch) = ephem_novatel2gavlab(ephem_novatel(:,ch));
  
  for t = 1:length(svpos_gpstimes)
    [ svpos{t,ch}, svvel{t,ch}, svclkcorr(t,ch) ] = ...
      calc_sv_pos_and_vel( ephem_novatel(:,ch), svpos_gpstimes(t), transit_time(t,ch) );
  end
  
end
clear eph_

userpos0_ecef = wgslla2xyz( gNovatel.zLat(1), gNovatel.zLong(1), gNovatel.zHeight(1) );

% Plot positions
for k = 1:nsv_have_ephem
  prnstr{k} = num2str(svs_have_ephem(k));
end
plot_svpos(svpos,userpos0_ecef,prnstr);
clear k prnstr


%% Do Acquisition
% This needs to output a struct that contains the same data as `acq` in the file
% `../lab4/part2a_narrow_ack.m`

acq_ugh = load('run_north_high_res_acq.mat');
acq = acq_ugh.run_north_high_res_acq;
clear acq_ugh
  


%% Do tracking with Akos' version to get IP

trackRes_akos = akos_tracking(nordnav_fid, acq, stop_time);
save('run_north_high_res_trackRes_akos.mat','trackRes_akos');
save workspace_after_akos_track


%% Do time syncing

clear all; close all; clc
load workspace_after_akos_track

% add dummy sv
trackRes_akos(5) = trackRes_akos(4);

gpsseconds_week = 60*60*24*7;

% number of milliseconds of data (CA code periods)
n_code_per = length(trackRes_akos(1).IP);
for k = 1:length(trackRes_akos)
  if length(trackRes_akos(k).IP) ~= n_code_per
    error('Nonequal IP lengths');
  end
end

[TLM_starts, TOW, TOW_empty, ch_status] = getTOWatIPIdx(trackRes_akos, acq, n_code_per);

% % Find SV positions at each GPS time reported.
% Each SV should have the same set of GPS times, just at different indices
% within the IP.
% get data from Novatel

% Calculate TOF at each GPS time 
transit_time_nordnav = zeros(length(TOW),4);
% calculate the transit time using Novatel's PSRs
for ch = 1:4
  transit_time_nordnav(:,ch) = interp1(svpos_gpstimes, transit_time(:,acq.svs(ch)), TOW(:,ch), 'cubic');
end

% get range of GPS times for which to calculate the SV Positions

  
% sv positions at those GPS times
%   - rows are each time, cols are each sv
svpos_nordnav = cell(length(TOW), 4); 
svvel_nordnav = svpos_nordnav;
% SV clock corrections, dimension correspond to `svpos`
svclkcorr_nordnav = zeros(size(svpos_nordnav));


for ch = 1:4
    
  for t = 1:length(TOW)
    [ svpos_nordnav{t,ch}, svvel_nordnav{t,ch}, svclkcorr_nordnav(t,ch) ] = ...
      calc_sv_pos_and_vel( ephem_novatel(:,ch), TOW(t,ch), transit_time_nordnav(t,ch) );
  end
  
end
clear eph_
% % Plot positions
% for k = 1:4
%   prnstr{k} = num2str(acq.svs(k));
% end
% plot_svpos(svpos_nordnav,userpos0_ecef,prnstr);


% % subtract TOF from each known
% TOW = TOW + transit_time_nordnav - svclkcorr_nordnav;

% interpolate Xbow values at its GPS times to find Xbow values at IP GPS times
XBow_accY=zeros(size(TOW));

for k = 1:4
  XBow_accY(:,k)  = spline(gXbow440.gps_time/1000, gXbow440.zAccelY, TOW(:,k));
end

save workspace_after_time_sync




%% End Matters

% save data with the filter order and the run name in the file
% save(['output_run_' run '_filter_' num2str(filter_order) '.mat']);
toc



