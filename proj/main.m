%% main.m Primary File for Final Project
% 
% 
clear all; close all; clc
tic

%% Settings / Constants

% run = 'south';
run = 'north';

% filter_order = 1;
filter_order = 2;
% filter_order = 3;

data_dir = ['..' filesep 'data' filesep 'final_proj_data' filesep];

stop_time = 120; % how long to calculate in sec

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
  svpos_gpstimes(data_count) = gNovatel.zPsrL1{k}(34)/1000;
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


%% Do Acquisition
% This needs to output a struct that contains the same data as `acq` in the file
% `../lab4/part2a_narrow_ack.m`




%% Do tracking with Akos' version to get IP

% acq = [];
% trackRes = akos_tracking(nordnav_fid, acq, stop_time);


%% Do time syncing

% % number of milliseconds of data (CA code periods)
% n_code_per = length(trackRes(1).IP);
% for k = 1:length(trackRes)
%   if length(trackRes(k).IP) ~= n_code_per
%     error('Nonequal IP lengths');
%   end
% end
% 
% [TLM_starts, TOW, ch_status] = getTOWatIPIdx(trackRes, acq, n_code_per);

% % Find SV positions at each GPS time reported.
% Each SV should have the same set of GPS times, just at different indices
% within the IP.

% Calculate TOF at each GPS time 

% subtract TOF from each known

% interpolate to find GPS time at each IP data epoch

% interpolate Xbow values at its GPS times to find Xbow values at IP GPS times


%% Our Tracking comparison



%% End Matters

% save data with the filter order and the run name in the file
save(['output_run_' run '_filter_' num2str(filter_order) '.mat']);
toc



