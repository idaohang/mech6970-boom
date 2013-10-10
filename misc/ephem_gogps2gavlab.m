function [out, prn, t_oe] = ephem_gogps2gavlab(in)

% EPHEM_GOGPS2GAVLAB takes in an ephemeris navigation data matrix output by
% the goGPS function RINEX_get_nav and rearranges it into the format
% expected by the GAVLab function calc_sv_pos
% 
% INPUTS:
% 
%   in: ephemeris matrix as output by goGPS function RINEX_get_nav
%       Indices:
%           Eph(1,i)  = svprn;
%           Eph(2,i)  = af2;
%           Eph(3,i)  = M0;
%           Eph(4,i)  = roota;
%           Eph(5,i)  = deltan;
%           Eph(6,i)  = ecc;
%           Eph(7,i)  = omega;
%           Eph(8,i)  = cuc;
%           Eph(9,i)  = cus;
%           Eph(10,i) = crc;
%           Eph(11,i) = crs;
%           Eph(12,i) = i0;
%           Eph(13,i) = idot;
%           Eph(14,i) = cic;
%           Eph(15,i) = cis;
%           Eph(16,i) = Omega0;
%           Eph(17,i) = Omegadot;
%           Eph(18,i) = toe;
%           Eph(19,i) = af0;
%           Eph(20,i) = af1;
%           Eph(21,i) = toc;
%           Eph(22,i) = IODE;
%           Eph(23,i) = code_on_L2;
%           Eph(24,i) = weekno;
%           Eph(25,i) = L2flag;
%           Eph(26,i) = svaccur;
%           Eph(27,i) = svhealth;
%           Eph(28,i) = tgd;
%           Eph(29,i) = fit_int;
%           Eph(30,i) = (sys_index-1) + svprn; %satellite index (consistent with other observation arrays)
%           Eph(31,i) = int8(sys_id);
%           Eph(32,i) = weektow2time(weekno, toe, sys_id);
%           Eph(33,i) = weektow2time(weekno, toc, sys_id);
%   
% 
% OUTPUTS:
% 
%   out: ephemeris matrix as expected by GAVLab function calc_sv_pos
%   prn: SV # (1-32) corresponds to columns of out & in

insz = size(in);
out = zeros(21,insz(2));
prn = in(1,:);
t_oe = in(18,:);

rearrange = [...
  28
  21
  2
  20
  19
  10
  11
  8
  9
  14
  15
  5
  3
  6
  4
  18
  16
  12
  7
  17
  13
];

out = in(rearrange,:);

end