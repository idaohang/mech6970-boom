function lon = wrapTo180(lon)
%wrapTo180 Wrap angle in degrees to [-180 180]
%
%   lonWrapped = wrapTo180(LON) wraps angles in LON, in degrees, to the
%   interval [-180 180] such that 180 maps to 180 and -180 maps to -180.
%   (In general, odd, positive multiples of 180 map to 180 and odd,
%   negative multiples of 180 map to -180.)
%
%   See also wrapTo360, wrapTo2Pi, wrapToPi.

% Copyright 2007-2008 The MathWorks, Inc.
% $Revision: 1.1.6.1 $  $Date: 2008/12/22 23:50:52 $

q = (lon < -180) | (180 < lon);
lon(q) = wrapTo360(lon(q) + 180) - 180;
