%% Data Input
load(['..' filesep 'data' filesep 'Novatel_Data_ephemfixed.mat'])
clear time 

psrL1r0 = gNovatel0.zPsrL1_gNovatel0;
psrL1r1 = gNovatel1.zPsrL1_gNovatel1;
adrL1r0 = gNovatel0.zAdrL1_gNovatel0;
adrL1r1 = gNovatel1.zAdrL1_gNovatel1;
% adrL2r0 = gNovatel0.zAdrL2_gNovatel0;

ndat_r0 = length(psrL1r0); % number of measurement epochs
ndat_r1 = length(psrL1r1);
valid_dat_r0 = [];
valid_dat_r1 = [];

% Find the GPS milliseconds into GPS week for which we have raw measurements
time_r0 = zeros(1,length(psrL1r0));
time_r1 = zeros(1,length(psrL1r1));
for i = 1:ndat_r0
  if ~(isempty(psrL1r0{i}) || isempty(adrL1r0{i}))
    time_r0(i) = psrL1r0{i}(34);
    valid_dat_r0(end+1) = i;
    %!! assume that all ADR and PSR came in at the same time
  end
end
for i = 1:ndat_r1
  if ~(isempty(psrL1r1{i}) || isempty(adrL1r1{i}))
    time_r1(i) = psrL1r1{i}(34);
    valid_dat_r1(end+1) = i;
  end
end
% fixing measurements - remove invalid data epochs
time_r0 = time_r0(valid_dat_r0);
time_r1 = time_r1(valid_dat_r1);
psrL1r0_ = psrL1r0(valid_dat_r0);
psrL1r1_ = psrL1r1(valid_dat_r1);
adrL1r0_ = adrL1r0(valid_dat_r0);
adrL1r1_ = adrL1r1(valid_dat_r1);
ndat_r0 = length(time_r0);
ndat_r1 = length(time_r1);

% SV's for which ephemeris data exists
prns_r0 = [1 2 4 8 9 10 12 17 20 24 28 32]; % hard-coded
prns_r1 = prns_r0; % need them to match
nsv_r0 = length(prns_r0);
nsv_r1 = length(prns_r1);

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


% Turn PSR/ADR data into a matrix
psrL1r0 = zeros(nsv_r0,ndat_r0);
psrL1r1 = zeros(nsv_r1,ndat_r1);
adrL1r0 = zeros(nsv_r0,ndat_r0);
adrL1r1 = zeros(nsv_r1,ndat_r1);
for k = 1:ndat_r0
  psrL1r0(:,k) = psrL1r0_{k}(prns_r0+1);
  adrL1r0(:,k) = adrL1r0_{k}(prns_r0+1);
end
for k = 1:ndat_r1
  psrL1r1(:,k) = psrL1r1_{k}(prns_r0+1);
  adrL1r1(:,k) = adrL1r1_{k}(prns_r0+1);
end

% if there's no data, remove the sv
have_dat_r0 = find(any(psrL1r0,2));
have_dat_r1 = find(any(psrL1r1,2));
prns_r0 = prns_r0(have_dat_r0);
prns_r1 = prns_r1(have_dat_r1);
nsv_r0 = length(prns_r0);
nsv_r1 = length(prns_r1);
psrL1r0 = psrL1r0(have_dat_r0,:);
psrL1r1 = psrL1r1(have_dat_r1,:);
adrL1r0 = adrL1r0(have_dat_r0,:);
adrL1r1 = adrL1r1(have_dat_r1,:);
ephem_mat_novatel = ephem_mat_novatel(:,have_dat_r0);

% put time into seconds from milliseconds
time_r0 = time_r0/1000;
time_r1 = time_r1/1000;

% remove the extra measurement on the beginning or r0
time_r0 = time_r0(2:end);
adrL1r0 = adrL1r0(:,2:end);
psrL1r0 = psrL1r0(:,2:end);
ndat_r0 = ndat_r0 - 1;

% reconcile common stuff
time = time_r0;
nsv = nsv_r0; % same number for both.
prns = prns_r0;
ndat = ndat_r0;


% !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
% Bad coding ahead:
% Simply leave out SV 12 since it comes in part of the way through the data set.
psrL1r0 = psrL1r0([1:3,5:end],:);
psrL1r1 = psrL1r1([1:3,5:end],:);
adrL1r0 = adrL1r0([1:3,5:end],:);
adrL1r1 = adrL1r1([1:3,5:end],:);
prns = prns([1:3,5:end]);
nsv = nsv-1;
% !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


% PRN Labels (SV #)
for k = 1:length(prns)
  prns_label{k} = num2str(prns(k));
end


% clean up
clear valid_dat_r0 valid_dat_r1 adrL1_ adrL2_ psrL1_ have_dat_r0 have_dat_r1 
clear adrL1r0_ adrL1r1_ psrL1r0_ psrL1r1_ i k time_r0 time_r1 prns_r1 prns_r0
clear ndat_r0 ndat_r1 nsv_r0 nsv_r1
