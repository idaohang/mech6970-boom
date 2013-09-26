function [satPositions,satClkCorr]=calc_sv_pos(ephem_data,transmitTime,transitTime)

% CALC_SV_POS calculates satellite positions in ECEF coordinates
% 
%   INPUTS:
% 
%     ephem_data: the large matrix of ephemeris parameters
%       dimensions: 21xn , where n is the number of sv's (currently only
%           n=1 is supported)
%     transmitTime:
%     transitTime:
% 
%   OUTPUTS:
% 
%     satPositions:
%     satClkCorr:
% 

gpsPi = 3.1415926535898;  % Pi used in the GPS coordinate system

% The columns should correspond to these (if made multi-D)
T_GD=ephem_data(1);
t_oc=ephem_data(2);
a_f2=ephem_data(3);
a_f1=ephem_data(4);
a_f0=ephem_data(5);

C_rc=ephem_data(6);
C_rs=ephem_data(7);
C_uc=ephem_data(8);
C_us=ephem_data(9);
C_ic=ephem_data(10);
C_is=ephem_data(11);

Delta_n=ephem_data(12);
M_0=ephem_data(13);
e=ephem_data(14);
sqrt_A=ephem_data(15);
t_oe=ephem_data(16);
Omega_0=ephem_data(17);
i_0=ephem_data(18);
omega=ephem_data(19);
dot_Omega=ephem_data(20);
Idot=ephem_data(21);

%--- Constants for satellite position calculation -------------------------
Omegae_dot = 7.2921151467e-5;  % Earth rotation rate, [rad/s]
GM = 3.986005e14;      % Earth's universal [m^3/s^2]
F = -4.442807633e-10; % Constant, [sec/(meter)^(1/2)]

% Process each satellite =================================================
satNr = 1 ;

%--- Find time difference ---------------------------------------------
dt = check_t(transmitTime - t_oc);

%--- Calculate clock correction ---------------------------------------
satClkCorr(satNr) = (a_f2 * dt + a_f1) * dt + a_f0 - T_GD;
time = transmitTime - satClkCorr(satNr);

% Find satellite's position ----------------------------------------------
%Restore semi-major axis
a   = sqrt_A*sqrt_A;
%Time correction
tk  = check_t(time - t_oe);
%Initial mean motion
n0  = sqrt(GM / a^3);
%Mean motion
n   = n0 + Delta_n;
%Mean anomaly
M   = M_0 + n * tk;
%Reduce mean anomaly to between 0 and 360 deg
M   = rem(M + 2*gpsPi, 2*gpsPi);
%Initial guess of eccentric anomaly
E   = M;

%--- Iteratively compute eccentric anomaly ----------------------------
for ii = 1:10
	E_old = E;
	E = M + e * sin(E);
	dE = rem(E - E_old, 2*gpsPi);
	if abs(dE) < 1.e-12
		break; % Necessary precision is reached, exit from the loop
	end
end

%Reduce eccentric anomaly to between 0 and 360 deg
E = rem(E + 2*gpsPi, 2*gpsPi);
%Compute relativistic correction term
dtr = F * e * sqrt_A * sin(E);
%Calculate the true anomaly
nu = atan2(sqrt(1 - e^2) * sin(E), cos(E)-e);
%Compute angle phi
phi = nu + omega;
%Reduce phi to between 0 and 360 deg
phi = rem(phi, 2*gpsPi);

%Correct argument of latitude
u = phi + C_uc * cos(2*phi) + C_us * sin(2*phi);
%Correct radius
r = a * (1 - e*cos(E)) + C_rc * cos(2*phi) + C_rs * sin(2*phi);
%Correct inclination
i = i_0 + Idot * tk + C_ic * cos(2*phi) + C_is * sin(2*phi);

%Compute the angle between the ascending node and the Greenwich meridian
Omega = Omega_0 + (dot_Omega-Omegae_dot)*tk - Omegae_dot*t_oe - Omegae_dot*transitTime;
%Reduce to between 0 and 360 deg
Omega = rem(Omega + 2*gpsPi, 2*gpsPi);

%--- Compute satellite coordinates ------------------------------------
satPositions(1, satNr) = cos(u)*r * cos(Omega) - sin(u)*r * cos(i)*sin(Omega);
satPositions(2, satNr) = cos(u)*r * sin(Omega) + sin(u)*r * cos(i)*cos(Omega);
satPositions(3, satNr) = sin(u)*r * sin(i);

% Include relativistic correction in clock correction -----------------
satClkCorr(satNr) = (a_f2 * dt + a_f1) * dt + a_f0 - T_GD + dtr;

end


%%%%%%%%%%%%%%%%%%%%%%%%%%% Auxiliary Function %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function corrTime = check_t(time)
% CHECK_T accounting for beginning or end of week crossover.
%
% corrTime = check_t(time);
%
%   Inputs:
%       time        - time in seconds
%
%   Outputs:
%       corrTime    - corrected time (seconds)
% Kai Borre 04-01-96
% Copyright (c) by Kai Borre
%
% CVS record:
% $Id: check_t.m,v 1.1.1.1.2.4 2006/08/22 13:45:59 dpl Exp $
%==================================================================
half_week = 302400;     % seconds
corrTime = time;

if time > half_week
	corrTime = time - 2*half_week;
elseif time < -half_week
	corrTime = time + 2*half_week;
end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% end check_t.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


