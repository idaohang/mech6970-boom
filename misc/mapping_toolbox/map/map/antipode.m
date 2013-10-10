function [lat,lon] = antipode(lat,lon,units)
%ANTIPODE Point on opposite side of globe
%
%   [NEWLAT, NEWLON] = ANTIPODE(LAT, LON) returns the geographic
%   coordinates of the points exactly opposite on the globe from the
%   input points given by LAT and LON.  All angles are in degrees.
%
%   [NEWLAT, NEWLON] = ANTIPODE(LAT, LON, ANGLEUNITS) specifies the
%   input and output units with the string ANGLEUNITS.  ANGLEUNITS can
%   be either 'degrees' or 'radians'.  It may be abbreviated and is
%   case-insensitive.

% Copyright 1996-2007 The MathWorks, Inc.
% $Revision: 1.9.4.5 $  $Date: 2007/05/10 13:47:24 $

% Note:  In order to map a longitude of zero to 180 degrees
% (or to pi, when working in radians), we use formulas like:
%
%   lon = 180 - mod(-lon, 360);
%
% rather than the more obvious:
%
%   lon = mod(lon, 360) - 180;

lat = -lat;
if nargin < 3
    lon = 180 - mod(-lon, 360);
else
    if strncmpi(units,'radians',numel(units))
        lon = pi - mod(-lon, 2*pi);
    elseif strncmpi(units,'degrees',numel(units))
        lon = 180 - mod(-lon, 360);
    else
        checkAngleUnits(units)
    end
end
