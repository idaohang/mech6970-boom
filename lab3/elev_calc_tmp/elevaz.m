function [Elev,Az] = elevaz(R_sc)
R_sc = Rsc/1000; % put 
% R_sc    = [-2000;4500;-4500];
H       = 0.42;
lat     = 40.5;
lst     = 90.5;

Re      = 6378.137;     % Equatorial Earh's radius [km]
Rp      = 6356.7523;    % Polar Earh's radius [km]
f       = (Re - Rp)/Re; % Oblateness or flattening

C1   = (Re/(1 - (2*f - f^2)*sind(lat)^2)^0.5 + H)*cosd(lat);
C2   = (Re*(1 - f)^2/(1 - (2*f - f^2)*sind(lat)^2)^0.5 + H)*sind(lat);
% Position vector of the observer,GEF
R_ob = [C1*cosd(lst); C1*sind(lst);C2];
% Position vector of the spacecraft relative to the observer
R_rel = R_sc - R_ob;

% GE_TH is direction cosine matrix to transform position vector components
% from geocentric equatorial frame into the topocentric horizon fream

GE_TH = [-sind(lst)          cosd(lst)              0;
       -sind(lat)*cosd(lst) -sind(lat)*sind(lst)  cosd(lat);
        cosd(lat)*cosd(lst)  cosd(lat)*sind(lst)   sind(lat)
   ];
R_rel_TH = GE_TH*R_rel;
rv = R_rel_TH/norm(R_rel_TH);
Elev = asin(rv(3))*180/pi;      % Elevation angle
Az  =atan2(rv(1),rv(2))*180/pi; % Azimuth angle
% fprintf('Elevation angle =  %4.2f [deg] \n',Elev);
% fprintf('Azimuth angle   =  %4.2f [deg] \n',Az);