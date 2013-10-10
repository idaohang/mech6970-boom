function sm = nm2sm(nm)
%NM2SM Convert distance from nautical to statute miles
%
%  sm = NM2SM(nm) converts distances from nautical miles to statute miles.
%
%  See also SM2NM, NM2DEG, NM2RAD, NM2KM.

% Copyright 1996-2007 The MathWorks, Inc.
% $Revision: 1.9.4.4 $  $Date: 2007/03/27 19:12:26 $

% Exact conversion factor
% 1 nm = 1852 meters, 1200/3937 meters = 1 statute foot,
% 5280 statute feet = 1 statute mile
cf = 1*1852/(1200/3937)/5280;
sm = nm * cf;
