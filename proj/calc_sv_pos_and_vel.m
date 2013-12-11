function [SV_Pos,SV_Vel,SV_Clk_corr]=calc_sv_pos_and_vel(ephem_data,transmitTime,transitTime)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                Constants
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

T_GD=ephem_data(25);
t_oc=ephem_data(24);
a_f2=ephem_data(28);
a_f1=ephem_data(27);
a_f0=ephem_data(26);
C_rc=ephem_data(15);
C_rs=ephem_data(16);
C_uc=ephem_data(13);
C_us=ephem_data(14);
C_ic=ephem_data(17);
C_is=ephem_data(18);
Delta_n=ephem_data(9);
M_0=ephem_data(10);
e=ephem_data(11);
Sqrt_semiMajAx=sqrt(ephem_data(8));  
t_oe=ephem_data(7);
Omega_0=ephem_data(21);
i_0=ephem_data(19);
omega=ephem_data(12);
dot_Omega=ephem_data(22);
Idot=ephem_data(20);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                Constants
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
gpsPi = 3.1415926535898;  % Pi used in the GPS coordinate system
Omegae_dot = 7.2921151467e-5;  % Earth rotation rate (rad/s)
GM = 3.986005e14;      %  (m^3/s^2)
F = -4.442807633e-10; % Relativistic Constant, (sec/m^(1/2))


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                Fix Time
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%---Find time differance by removing sv clock correction term---% 
dt = check_t(transmitTime - t_oc); 

%---Calculate clock correction term---%
SV_Clk_corr= (a_f2 * dt + a_f1) * dt + a_f0 - T_GD;

tk = check_t(transmitTime - SV_Clk_corr - t_oe);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                           Find SV Position
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Restore semi-major axis
a=Sqrt_semiMajAx*Sqrt_semiMajAx;

%Initial mean motion
n0=sqrt(GM/a^3);

%Mean motion
n=n0 + Delta_n;

%Mean anomaly
M=M_0 + n * tk;

%Reduce mean anomaly to between 0 and 360 deg
M=rem(M + 2*gpsPi,2*gpsPi);

%Initial guess of eccentric anomaly
E   = M;

%Find Eccentric Anomaly
for i = 1:10
	E_old = E;
	E=M+e*sin(E);
	delE = rem(E-E_old,2*gpsPi);

	if abs(delE)<1.e-12
		break
	end
end

%Reduce eccentric anomaly to between 0 and 360 deg
E=rem(E + 2*gpsPi, 2*gpsPi);

%Compute relativistic correction term
dtr=F*e*Sqrt_semiMajAx*sin(E);

%Calculate the true anomaly
nu=atan2(sqrt(1 - e^2)*sin(E),cos(E)-e);

%Compute angle phi
phi=nu+omega;

%Reduce phi to between 0 and 360 deg
phi=rem(phi,2*gpsPi);

%Correct argument of latitude
u=phi+C_uc*cos(2*phi)+C_us*sin(2*phi);

%Correct radius
r=a*(1-e*cos(E))+C_rc*cos(2*phi)+C_rs*sin(2*phi);

%Correct inclination
i=i_0+Idot*tk+C_ic*cos(2*phi)+C_is*sin(2*phi);

%Compute the angle between the ascending node and the Greenwich meridian
Omega=Omega_0+(dot_Omega-Omegae_dot)*tk-Omegae_dot*t_oe-Omegae_dot*transitTime;

%Reduce to between 0 and 360 deg
Omega=rem(Omega+2*gpsPi,2*gpsPi);

X=r*cos(u);
Y=r*sin(u);

x=X*cos(Omega)-Y*cos(i)*sin(Omega);
y=X*sin(Omega)+Y*cos(i)*cos(Omega);
z=Y*sin(i);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                           Find SV Position
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Edot=(n0+Delta_n)/(1-e*cos(E));

phidot=(sqrt(1-e^2)/(1-e*cos(E)))*Edot;

udot=(1+2*C_us*cos(2*phi)-2*C_uc*sin(2*phi))*phidot;

rdot=2*(C_rs*cos(2*phi)-C_rc*sin(2*phi))*phidot + a*e*sin(E)*Edot;

idot=2*(C_is*cos(2*phi)-C_ic*sin(2*phi))*phidot + Idot;

Xdot=rdot*cos(u) - r*sin(u)*udot;

Ydot=rdot*sin(u) + r*cos(u)*udot;

Omegadot=dot_Omega-Omegae_dot;

xdot=Xdot*cos(Omega)-Ydot*cos(i)*sin(Omega)+Y*sin(i)*sin(Omega)*idot-y*Omegadot;

ydot=Xdot*sin(Omega)+Ydot*cos(i)*cos(Omega)-Y*sin(i)*cos(Omega)*idot+x*Omegadot;

zdot=Ydot*sin(i)+Y*cos(i)*idot;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                           Creating Outputs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

SV_Pos=[x,y,z];
SV_Vel=[xdot,ydot,zdot];
SV_Clk_corr = (a_f2 * dt + a_f1) * dt + a_f0 - T_GD + dtr;

end





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                          Some Other Function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function corrTime = check_t(time)
%CHECK_T accounting for beginning or end of week crossover.
%
%corrTime = check_t(time);
%
%   Inputs:
%       time        - time in seconds
%
%   Outputs:
%       corrTime    - corrected time (seconds)

%Kai Borre 04-01-96
%Copyright (c) by Kai Borre
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