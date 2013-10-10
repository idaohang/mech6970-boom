function [out,prns,t_ow] = ephem_novatel2gavlab(in)

% DESCRIPTION: given a ephemeris data as parsed out of a novatel alog,
% return something in the order expected by calc_sv_pos
% 
% INPUTS:
% 
%   in: 32x30 matrix of ephemeris data in the order output by novatel zEphemNN
%     rows correspond to SV, columns correspond to data fields
%       Indices:
%         1:  time of week
%         2:  health
%         3:  issue of ephem data 1
%         4:  issue of ephem data 2
%         5:  week
%         6:  z week
%         7:  toe
%         8:  A, semimajor axis
%         9:  delta N, mean motion diff, rad/s
%         10: M_0, mean anomaly of reference time
%         11: eccentricity
%         12: omega - argument of perigee
%         13: cuc
%         14: cus
%         15: crc
%         16: crs
%         17: cic
%         18: cis
%         19: I_0, inclination angle
%         20: I_dot, rate of inclination angle
%         21: omega_0, right ascension
%         22: omega_dot, rate of right ascension
%         23: iodc, issue of data clock
%         24: toc, sv clock correction term
%         25: tgd, group delay difference
%         26: af0
%         27: af1
%         28: af2
%         29: cmot, corrected mean motion, rad/s
%         30: ura, user range accuracy variance,  m^2
% 
% 
% OUTPUTS:
% 
%   out: see input `in` of gps.calc_sv_pos
%     21xn matrix, columns are each data field, rows are each SV
%     Indices:
%         T_GD=ephem_data(1);
%         t_oc=ephem_data(2);
%         a_f2=ephem_data(3);
%         a_f1=ephem_data(4);
%         a_f0=ephem_data(5);
%         C_rc=ephem_data(6);
%         C_rs=ephem_data(7);
%         C_uc=ephem_data(8);
%         C_us=ephem_data(9);
%         C_ic=ephem_data(10);
%         C_is=ephem_data(11);
%         Delta_n=ephem_data(12);
%         M_0=ephem_data(13);
%         e=ephem_data(14);
%         sqrt_A=ephem_data(15);
%         t_oe=ephem_data(16);
%         Omega_0=ephem_data(17);
%         i_0=ephem_data(18);
%         omega=ephem_data(19);
%         dot_Omega=ephem_data(20);
%         Idot=ephem_data(21);
% 
%   prns: which prns the columns correspond to
%   t_ow: time of week for ephemeris

have_data = in(:,1)~=0; % indices correspond to SV# for which data exists
[prns,~] = find(have_data);

rearrange = [...
  25
  24
  28
  27
  26
  15
  16
  13
  14
  17
  18
  9
  10
  11
  8 % semimajor axis, not the sqrt
  7
  21
  19
  12
  22
  20
];

% out = zeros(21,length(in(have_data,1))); 
out = in(prns,rearrange);
% make A -> sqrt(A)
out(:,15) = sqrt(out(:,15));
% output time of ephemeris for each SV's data
t_ow = in(prns,1);

end



