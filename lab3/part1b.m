%% MECH 6970 Lab 3, Part I, b)
% 
% Robert Cofield
% 
% created 2013-10-09
% Using Receiver 0
%   had to copy some ephemeris data from Receiver 1 over to Receiver 0 due to
%   parsing errors.
% 

%% Data Input
clear all; close all; clc
fprintf('\nPart 1 - b)\n')
load(['..' filesep 'data' filesep 'Novatel_Data_ephemfixed.mat'])
clear time gNovatel1

psrL1 = gNovatel0.zPsrL1_gNovatel0;
adrL1 = gNovatel0.zAdrL1_gNovatel0;
adrL2 = gNovatel0.zAdrL2_gNovatel0;

ndat = length(psrL1); % number of measurement epochs
valid_dat = [];

% Find the GPS milliseconds into GPS week for which we have raw measurements
time = zeros(1,length(psrL1));
for i = 1:ndat
  if isempty(psrL1{i}) || isempty(adrL2{i}) || isempty(adrL2{i})
    continue
  end
  time(i) = psrL1{i}(34);
  valid_dat(end+1) = i;
  %!! assume that all ADR and PSR came in at the same time
end
% fixing measurements - remove invalid data epochs
time = time(valid_dat);
psrL1_ = psrL1(valid_dat);
adrL1_ = adrL1(valid_dat);
adrL2_ = adrL2(valid_dat);
ndat = length(time);

% SV's for which ephemeris data exists
prns = [1 2 4 8 9 10 12 17 20 24 28 32];
nsv = length(prns);

% Turn PSR/ADR data into a matrix
psrL1 = zeros(nsv,ndat);
adrL1 = zeros(nsv,ndat);
adrL2 = zeros(nsv,ndat);
for k = 1:ndat
  psrL1(:,k) = psrL1_{k}(prns+1);
  adrL1(:,k) = adrL1_{k}(prns+1);
  adrL2(:,k) = adrL2_{k}(prns+1);
end

% if there's no data, remove the sv
have_dat = find(any(psrL1,2));
prns = prns(have_dat);
nsv = length(prns);
psrL1 = psrL1(have_dat,:);
adrL1 = adrL1(have_dat,:);
adrL2 = adrL2(have_dat,:);

clear valid_dat adrL1_ adrL2_ psrL1_ have_dat

%% Parameters

% Smoothing window (samples @ 1 Hz)
M = 100;

% Estimated transit time from SV to user
c = 299792458;
transit_time_est = 20e6/c;

% ECEF estimate of the user position for unit vectors (using Google Earth)
lla_user_est = [dms2degrees(32,35,26.1), -dms2degrees(85,29,20.61), 205]; % lat lon alt
user_pos_est = 


%% Calculate SV Positions from Ephemeris


% ephem_mat: each column corresponds to an sv
ephem_mat_novatel = [...
  gNovatel0.zEphem1_gNovatel0{end},...
  gNovatel0.zEphem2_gNovatel0{end},...
  gNovatel0.zEphem4_gNovatel0{end},...
  gNovatel0.zEphem8_gNovatel0{end},...
  gNovatel0.zEphem9_gNovatel0{end},...
  gNovatel0.zEphem10_gNovatel0{end},...
  gNovatel0.zEphem12_gNovatel0{end},...
  gNovatel0.zEphem17_gNovatel0{end},...
  gNovatel0.zEphem20_gNovatel0{end},...
  gNovatel0.zEphem24_gNovatel0{end},...
  gNovatel0.zEphem28_gNovatel0{end},...
  gNovatel0.zEphem32_gNovatel0{end}...
];

ephem_mat = zeros(21,nsv);
ephem_time = zeros(1,nsv); % gps seconds into week at which ephem was transmitted
sv_clkcorr = zeros(nsv,ndat);
svpos = zeros(3,nsv,ndat);

for i = 1:nsv
  % get the ephemeris matrix into the correct format
  [ephem_mat(:,i), ephem_time(i)] = ephem_novatel2gavlab(ephem_mat_novatel(:,i));
  % Find the SV Position at each of the measurement epochs we have
  for k = 1:ndat
    [svpos(:,i,k), sv_clkcorr(i,k)] = calc_sv_pos(ephem_mat(:,i), time(k), transit_time_est);
  end
end

% Look at SV position LLA just to check the positions


clear ephem_mat_novatel i k ephem_mat

%% Single Frequency Carrier Smoothing
% 

% Carrier smoothed range estimates
r1 = zeros(nsv,ndat);
r1(:,1) = psrL1(:,1); % start by copying

for k = 2:ndat
  for i = 1:nsv
  
    % make sure we have data for this sv
    % assume that no PSR data means no ADR data, vice versa
    if ~psrL1(i,k-1) % haven't had data before this
      if ~psrL1(i,k) % still don't have dat
        continue
      else % sv just came into view -> copy to begin (this happens to SV 12)
        r1(i,k) = psrL1(i,k);
        continue
      end
    end
    
    r1(i,k) = psrL1(i,k)/M + (M-1)/M*( r1(i,k-1) + adrL1(i,k) - adrL1(i,k-1));

  end
end

% Least Squares Estimation



%% Dual Frequency Carrier Smoothing
% 






















