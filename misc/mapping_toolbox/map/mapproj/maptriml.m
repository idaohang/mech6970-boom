function [lat, lon] = maptriml(lat, lon, latlim, lonlim)
%MAPTRIML  Trim lines to latitude-longitude quadrangle
%
%   [latTrimmed, lonTrimmed] = MAPTRIML(LAT, LON, LATLIM, LONLIM) trims
%   a line with vertices specified by vectors LAT and LON to the
%   quadrangle specified by LATLIM and LONLIM.  LATLIM is a vector of the
%   form [southern-limit northern-limit] and LONLIM is a vector of the
%   form [western-limit eastern-limit].  All angles are in units of
%   degrees.  Outputs latTrimmed and lonTrimmed are column vectors
%   regardless of the shape of inputs LAT and LON.
%
%   See also MAPTRIMP, MAPTRIMS.

% Copyright 1996-2007 The MathWorks, Inc.
% $Revision: 1.1.6.1 $  $Date: 2007/11/26 20:35:57 $

checklatlon(lat, lon, mfilename, 'LAT', 'LON', 1, 2)
checkgeoquad(latlim, lonlim, mfilename, 'LATLIM', 'LONLIM', 3, 4)
[lat, lon, latlim, lonlim] = toRadians('degrees', lat, lon, latlim, lonlim);
[lat, lon] = trimPolylineToQuadrangle(lat, lon, latlim, lonlim);
[lat, lon] = toDegrees('radians', lat(:), lon(:));
