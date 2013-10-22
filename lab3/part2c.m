%% MECH 6970 Lab 3, Part 2, c)
% 
% Robert Cofield
% 
% Note that for this you need the lambda source code in the `misc` folder added
% to your path
clear all; close all; clc
tic
part2c_load_data
try
  matlabpool(3) % comment this out if you don't have the parallel toolbox
end


%% Parameters

c = 299792458; % Speed of light, m/s
fL1 = 1575.42e6; % L1 frequency, Hz
fL2 = 1227.60e6; % L2 frequency, Hz
wavelengthL1 = c/fL1;
wavelengthL2 = c/fL2;

lla_user_est = [dms2deg([32,35,26.1]), -dms2deg([85,29,20.61]), 205]; % lat lon alt


%% Calculate SV Positions from Ephemeris

ecef_user_est = wgslla2xyz(lla_user_est(1),lla_user_est(2),lla_user_est(3));

ephem_mat = zeros(21,nsv);
ephem_time = zeros(1,nsv); % gps seconds into subplot(week at which ephem was transmitted
sv_clkcorr = zeros(nsv,ndat);
sv_clkcorr_psr = zeros(size(sv_clkcorr)); % range correction corresponding to clock correction
svpos = zeros(3,nsv,ndat);
% svpos_psr = zeros(nsv,ndat); % pseudorange using estimated user pos and svpos - error checking
svelev = zeros(nsv,ndat);
psrL1r0corr = zeros(size(psrL1r0));
psrL1r1corr = zeros(size(psrL1r1));

transit_time_est = 20e6/c; % seconds

for i = 1:nsv
  % get the ephemeris matrix into the correct format
  [ephem_mat(:,i), ephem_time(i)] = ephem_novatel2gavlab(ephem_mat_novatel(:,i));
  % Find the SV Position at each of the measurement epochs we have
  ephem_mat_ = ephem_mat(:,i);
    
  % find elevation of each satellite


    parfor k = 1:ndat
    tx_time = time(k)%;-transit_time_est;
    [svpos(:,i,k), sv_clkcorr(i,k)] = calc_sv_pos(ephem_mat_, tx_time, transit_time_est);
    sv_clkcorr_psr(i,k) = sv_clkcorr(i,k)*c;
    psrL1r0corr(i,k) = psrL1r0(i,k) + sv_clkcorr_psr(i,k);
    psrL1r1corr(i,k) = psrL1r1(i,k) + sv_clkcorr_psr(i,k);
    
%     % check position using approx psr
%     if psrL1r0(i,k)
%       svpos_psr(i,k) = norm( svpos(:,i,k)-ecef_user_est ,2);
%     end
    
  end
end

% find the elevation of each satellite at each epoch
parfor k = 1:ndat
  [~,svelev(:,k)] = calc_azel(ecef_user_est',reshape(svpos(:,:,k),nsv,3));
end


% % Find the elevation of each at the beginning
% ecef_user_est = wgslla2xyz(lla_user_est(1), lla_user_est(2), lla_user_est(3));
% svpos0 = reshape(svpos(:,:,end),8,3);
% [~,elevation] = calc_azel(ecef_user_est',svpos0);

% Look at SV position LLA just to check the positions
svpos_lla = zeros(size(svpos));
parfor i = 1:nsv
  for k = 1:ndat
    [lat, lon, alt] = coordutil.wgsxyz2lla(svpos(:,i,k),1000);
    svpos_lla(:,i,k) = [lat, lon, alt];
  end
end

% % output LLA sv pos initial to KML file
% kml_f_svpos0 = kml_file.createFolder('SV Initial Positions');
% for k = 1:8
%   kml_f_svpos0.point(svpos_lla(2,k,1),svpos_lla(1,k,1),svpos_lla(3,k,1));
% end

clear  i k svpos0 ephem_mat_novatel

%% Single Differencing

% turn adr into car
carL1r0 = -adrL1r0*wavelengthL1;
carL1r1 = -adrL1r1*wavelengthL1;

carL1r01sd = carL1r1 - carL1r0;
psrL1r01sd = psrL1r1 - psrL1r0;
Nsd01_float = (carL1r01sd-psrL1r01sd)/wavelengthL1;


% Ambiguity covariance
% Qsd01 = 

% [afixed, sqnorm] = LAMBDA(

%% Plot Single Differencing Results

figure;
  plot(Nsd01_float');
  title('Single Differenced Floating Ambiguity Estimates');
  ylabel('Est N_{ab}'); xlabel('sample #');
  legend(prns_label); grid on


%% Double Differencing

% Base PRN is the one highest in the sky
% base_prn_idx = 









%% End matters 

try
  matlabpool close
end
toc

