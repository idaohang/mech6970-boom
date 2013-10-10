function lon = westof(lon, meridian, units)
%WESTOF Wrap longitudes to values west of specified meridian
% 
%   WESTOF is obsolete and will be removed in a future release of
%   Mapping Toolbox.  Replace it with the following calls, which are
%   also more efficient:
%
%      westof(lon, meridian, 'degrees') ==>
%                     meridian - mod(meridian - lon, 360)
% 
%      westof(lon, meridian, 'radians') ==>
%                      meridian - mod(meridian - lon, 2*pi)
%
%   NEWLON = WESTOF(LON, MERIDIAN) wraps angles in LON to values in
%   the interval (MERIDIAN-360 MERIDIAN].  LON is a scalar longitude or
%   vector of longitude values.  All inputs and outputs are in degrees.
%
%   NEWLON = WESTOF(LON, MERIDIAN, ANGLEUNITS) specifies the input
%   and output units with the string ANGLEUNITS.  ANGLEUNITS can be
%   either 'degrees' or 'radians'.  It may be abbreviated and is
%   case-insensitive.  If ANGLEUNITS is 'radians', then the input is
%   wrapped to the interval (MERIDIAN-2*pi MERIDIAN].

% Copyright 1996-2010 The MathWorks, Inc.
% $Revision: 1.5.4.6 $  $Date: 2010/03/04 16:22:22 $

warning('map:westof:obsolete', ...
    ['Function %s is obsolete and will be removed in a future release of the toolbox. \n', ...
    'See the help for a more efficient alternative.'],...
    upper(mfilename))

if nargin < 3
    lon = meridian - mod(meridian - lon, 360);
else
    if strncmpi(units,'radians',numel(units))
        lon = meridian - mod(meridian - lon, 2*pi);
    elseif strncmpi(units,'degrees',numel(units))
        lon = meridian - mod(meridian - lon, 360);
    else
        checkAngleUnits(units)
    end
end
