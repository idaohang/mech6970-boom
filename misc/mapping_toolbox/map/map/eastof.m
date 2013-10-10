function lon = eastof(lon,meridian,units)
%EASTOF Wrap longitudes to values east of specified meridian
%
%   EASTOF is obsolete and will be removed in a future release of
%   Mapping Toolbox.  Replace it with the following calls, which are
%   also more efficient:
%
%      eastof(lon, meridian, 'degrees') ==>
%                      meridian + mod(lon - meridian, 360)
% 
%      eastof(lon, meridian, 'radians') ==>
%                      meridian + mod(lon - meridian, 2*pi)
%
%   NEWLON = EASTOF(LON, MERIDIAN) wraps angles in LON to values in
%   the interval [MERIDIAN MERIDIAN+360).  LON is a scalar longitude or
%   vector of longitude values.  All inputs and outputs are in degrees.
%
%   NEWLON = EASTOF(LON, MERIDIAN, ANGLEUNITS) specifies the input
%   and output units with the string ANGLEUNITS.  ANGLEUNITS can be
%   either 'degrees' or 'radians'.  It may be abbreviated and is
%   case-insensitive.  If ANGLEUNITS is 'radians', then the input is
%   wrapped to the interval [MERIDIAN MERIDIAN+2*pi).

% Copyright 1996-2010 The MathWorks, Inc.
% $Revision: 1.5.4.6 $  $Date: 2010/03/04 16:22:13 $

warning('map:eastof:obsolete', ...
    ['Function %s is obsolete and will be removed in a future release of the toolbox. \n', ...
    'See the help for a more efficient alternative.'],...
    upper(mfilename))

if nargin < 3
    lon = meridian + mod(lon - meridian, 360);
else
    if strncmpi(units,'radians',numel(units))
        lon = meridian + mod(lon - meridian, 2*pi);
    elseif strncmpi(units,'degrees',numel(units))
        lon = meridian + mod(lon - meridian, 360);
    else
        checkAngleUnits(units)
    end
end
