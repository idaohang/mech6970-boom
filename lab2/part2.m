%% MECH 6970 Lab 2, Part 2
% 
genutil.ccc
gps.constants

user_lla = [dms2dd(32,35,26.16), -dms2dd(85,29,21.21), 200];
user_ecef = coordutil.wgslla2xyz(user_lla(1),user_lla(2),user_lla(3));

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
svpos_lla = zeros(3,numsat);
svpos_ae = zeros(2,numsat);
for k = 1:length(prns)
  prn = prns(k);
  [pos,~] = gps.calc_sv_pos(ephem(k,:), tow(k), t_transit_est);
  svpos(:,k) = pos;
  
  [sv_lat,sv_lon,sv_alt] = coordutil.wgsxyz2lla(pos,100);
  
  % SV pos relative to user in ENU
  [sv_lat,sv_lon]=wraplatlong(sv_lat,sv_lon);

  if abs(sv_alt)<17000000
      continue;
  end
  
  svpos_lla(:,k) = [sv_lon,sv_lat,sv_alt]';

  fprintf('sv positions prn %i,\n\t LLA: %20.10f\t%20.10f\t%20.10f\n', k,sv_lat,sv_lon,sv_alt);
  
  dp_enu = coordutil.wgslla2enu(sv_lat,sv_lon,sv_alt, user_lla(1),user_lla(2),user_lla(3));
  [a,e,r] = cart2sph(dp_enu(1),dp_enu(2),dp_enu(3)); 
  svpos_ae(:,prn) = [a;e];
  fprintf('\tAzi, Ele: %7.2f deg, %7.2f deg\n', rad2deg(a),rad2deg(e))
end

%% Skyview

% svpos_ae(svpos_ae==0) = [];
% svpos_ae = [svpos_ae(1:9);svpos_ae(10:end)];

% figure;
% skyplot(svpos',prns,user_ecef,0,1);

figure;
polar_data = zeros(size(svpos_ae));
polar(0,0,'k'); hold on
view([90 -90]);
for k = 1:length(svpos_ae)
  if ~svpos_ae(1,k)
    continue
  end
  
  polar_data(1,k) = rad2deg(svpos_ae(1,k));
  polar_data(2,k) = svpos_ae(2,k) / (pi/2);
end
polar(polar_data','b.')

%% KML Google Earth

gek = kml('Lab2');
f = gek.createFolder('Receiver Positions');

for k = 1:length(gNovatel0.X.val)
  ecef = [gNovatel0.X.val(k) gNovatel0.Y.val(k) gNovatel0.Z.val(k)];
  [lat, lon, alt] = coordutil.wgsxyz2lla(ecef);
  f.point(lon,lat,alt);
end
gek.run()

%% ENU from Toomer's Corner

lla_ref = [dms2dd(32,36,23.52), -dms2dd(85,28,54.26), 219];
enu = zeros(3,length(gNovatel0.X.val));
for k = 1:length(gNovatel0.X.val)
  ecef = [gNovatel0.X.val(k) gNovatel0.Y.val(k) gNovatel0.Z.val(k)];
  enu(:,k) = coordutil.wgsxyz2enu(ecef', lla_ref(1),lla_ref(2),lla_ref(3));
end

figure;
  subplot(3,1,1)
    plot(gNovatel0.X.moos_time, enu(1,:))
    grid on
    ylabel('East (m)')
    xlabel('Database time (s)')
  subplot(3,1,2)
    plot(gNovatel0.X.moos_time, enu(2,:))
    grid on
    ylabel('North (m)')
    xlabel('Database time (s)')
  subplot(3,1,3)
    plot(gNovatel0.X.moos_time, enu(3,:))
    grid on
    ylabel('Up (m)')
    xlabel('Database time (s)')
  
open('./part2_n3_ks_enu.fig')
  
  
  









