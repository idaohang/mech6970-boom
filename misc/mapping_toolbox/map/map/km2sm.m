function sm = km2sm(km)
%KM2SM Convert distance from kilometers to statute miles
%
%  sm = KM2SM(km) converts distances from kilometers to statute miles.
%
%  See also SM2KM, KM2DEG, KM2RAD, KM2NM.

% Copyright 1996-2007 The MathWorks, Inc.
% $Revision: 1.9.4.4 $  $Date: 2007/03/27 19:12:15 $

% Exact conversion factor
% 1 kilometer = 1000 meters, 1200/3937 meters = 1 statue foot,
% 5280 statute feet = 1 statue mile
cf = 1*1000/(1200/3937)/5280;
sm = cf * km;
