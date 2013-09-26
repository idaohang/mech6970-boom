%% MECH 6970 Lab 2, Part 2
% 
genutil.ccc
gps.constants

%% Load data
% 
load(['..' filesep 'data' filesep 'Novatel_Data__pyparsed.mat'])
% pick Novatel0 to use for data
clear gNovatel1

% put ephem into a matrix
% use the latest data for each, if multiple ephem packets
ephem_novatel = zeros(32,30);
ephem_novatel([1,2,4,8,9,12,17,24,28,32],:) = [...
  gNovatel0.Ephem1.val(end,:)
  gNovatel0.Ephem2.val(end,:)
  gNovatel0.Ephem4.val(end,:)
  gNovatel0.Ephem8.val(end,:)
  gNovatel0.Ephem9.val(end,:)
  gNovatel0.Ephem12.val(end,:)
  gNovatel0.Ephem17.val(end,:)
  gNovatel0.Ephem24.val(end,:)
  gNovatel0.Ephem28.val(end,:)
  gNovatel0.Ephem32.val(end,:)
];
% convert columns to the format expected by calc_sv_pos
[ephem, prns, tow] = gps.ephem_novatel2gavlab(ephem_novatel);
numsat = length(prns);
%% Skyview
% 

% transit time estimation
range_est = 20e6;
t_transit_est = range_est/c;

% position estimation
svpos = zeros(3,numsat);
svpos_ae = zeros(2,numsat);
for k = 1:length(prns)
  prn = prns(k);
  [pos,clkcorr] = gps.calc_sv_pos(ephem(k,:), tow(k), t_transit_est);
  svpos(:,k) = pos;
  
  [sv_lat,sv_lon,sv_alt] = coordutil.wgsxyz2lla(pos);
  % SV pos relative to user in ENU
  [sv_lat,sv_lon]=wraplatlong(sv_lat,sv_lon);
  if abs(sv_alt)<18000000
      continue;
  end
  fprintf('sv positions, LLA: %20.10f\t%20.10f\t%20.10f\n',sv_lat,sv_lon,sv_alt);
  
  dp_enu = coordutil.wgslla2enu(sv_lat,sv_lon,sv_alt, user_lla(1),user_lla(2),user_lla(3));
  [a,e,r] = cart2sph(dp_enu(1),dp_enu(2),dp_enu(3)); 
  svpos_ae(:,prn,prn_data_cnt(prn)) = [a;e];

end