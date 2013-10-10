function lon = smoothlong(lon, units)
%SMOOTHLONG Remove discontinuities in longitude data
%
%   SMOOTHLONG is obsolete and has been replaced by unwrapMultipart,
%   which requires its input to be in radians.  When working in degrees,
%   use radtodeg(unwrapMultipart(degtorad(lon))).
%
%   NEWLON = SMOOTHLONG(LON) unwraps a row or column vector of
%   longitudes, azimuths, or phase angles.  Input and output are both in
%   degrees.
%
%   NEWLON = SMOOTHLONG(LON, ANGLEUNITS) works in the units defined by
%   the string ANGLEUNITS, which can be either 'degrees' or 'radians'.
%   ANGLEUNITS may be abbreviated and is case-insensitive.

% Copyright 1996-2009 The MathWorks, Inc.
% $Revision: 1.5.4.6 $  $Date: 2009/05/14 17:05:33 $

warning('map:smoothlong:obsolete', ...
    ['Function %s is obsolete and has been replaced by unwrapMultipart, \n', ...
    'which requires its input to be in radians. When working in degrees, \n',...
    'use radtodeg(unwrapMultipart(degtorad(lon)))'], ...
    upper(mfilename))

if nargin < 2
    lon = radtodeg(unwrapMultipart(degtorad(lon)));
else
    if strncmpi(units,'radians',numel(units))
        lon = unwrapMultipart(lon);
    elseif strncmpi(units,'degrees',numel(units))
        lon = radtodeg(unwrapMultipart(degtorad(lon)));
    else
        checkAngleUnits(units)
    end
end
