%% MECH 6970 Lab 2, Part 2
% 
genutil.ccc
gps.constants

% solar house field
user_lla = [dms2dd(32,35,26.16), dms2dd(85,29,21.21), 200];

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
%% Satellite positions
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
  
  [sv_lat,sv_lon,sv_alt] = coordutil.wgsxyz2lla(pos,100);
  % SV pos relative to user in ENU
  [sv_lat,sv_lon]=wraplatlong(sv_lat,sv_lon);
  if abs(sv_alt)<18000000
      continue;
  end
  fprintf('sv positions prn %i,\n\t LLA: %20.10f\t%20.10f\t%20.10f\n', k,sv_lat,sv_lon,sv_alt);
  
  dp_enu = coordutil.wgslla2enu(sv_lat,sv_lon,sv_alt, user_lla(1),user_lla(2),user_lla(3));
  [a,e,r] = cart2sph(dp_enu(1),dp_enu(2),dp_enu(3)); 
  svpos_ae(:,prn) = [a;e];
  fprintf('\tAzi, Ele: %7.2f deg, %7.2f deg\n', rad2deg(a),rad2deg(e))
end

%% Skyview

figure;
polar(0,0, 'k.')
view([90 -90]); % get correct rotation style
hold on;
for k = 1:32
  if ~svpos_ae(1,k) % no sv data
    continue
  end
  
  azi = svpos_ae(1,k);
  ele = svpos_ae(2,k);
  polar(rad2deg(azi),rad2deg(ele), 'b')
  
end

%% KML Google Earth

















