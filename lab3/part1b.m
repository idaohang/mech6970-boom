%% Lab 3, Part I, b)
% Robert Cofield
% created 2013-10-09
% 

%% Data Input
clear all; close all; clc
fprintf('\nPart 1 - d)\n')
load(['..' filesep 'data' filesep 'Novatel_Data.mat'])
clear time gNovatel1

time = gNovatel0.time;
psrL1 = gNovatel0.zPsrL1_gNovatel0;
adrL1 = gNovatel0.zAdrL1_gNovatel0;
adrL2 = gNovatel0.zAdrL2_gNovatel0;

%% Parameters

% Smoothing window (samples @ 1 Hz)
M = 100;

% Estimated transit time from SV to user
c = 299792458;
transit_time_est = 20e6/c;


%% Calculate SV Positions from Ephemeris

% SV's for which ephemeris data exists
prns = [4 8 9 10 12 17 20 24 28 32];
nsv = length(prns);

% ephem_mat: each column corresponds to an sv
ephem_mat_novatel = [...
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
ephem_time = zeros(1,nsv);
sv_clkcorr = zeros(1,nsv);
svpos = zeros(3,nsv);

for i = 1:nsv
  [ephem_mat(:,i), ephem_time(i)] = ephem_novatel2gavlab(ephem_mat_novatel(:,i));
  [svpos(:,i), sv_clkcorr(i)] = calc_sv_pos(ephem_mat(:,i),ephem_time(i),transit_time_est);
end


%% Single Frequency Carrier Smoothing
% 




%% Dual Frequency Carrier Smoothing
% 