function [lat, lon] = maptrimp(lat, lon, latlim, lonlim)
%MAPTRIMP  Trim polygons to latitude-longitude quadrangle
%
%   [latTrimmed, lonTrimmed] = MAPTRIMP(lat, lon, latlim, lonlim) trims
%   the polygons in lat and lon to the quadrangle specified by latlim
%   and lonlim.  Latlim and lonlim are two element vectors, defining the
%   latitude and longitude limits respectively.  lat and lon must be
%   vectors that represent valid polygons.
%
%   MAPTRIMP conditions the longitude limits such that:
%     (1) lonlim(2) always exceeds lonlim(1)
%     (2) lonlim(2) never exceeds lonlim(1) by more than 360
%     (3) lonlim(1) < 180 or lonlim(2) > -180
%     (4) ensure that if the quadrangle spans the Greenwich meridian,
%         then that meridian appears at longitude == 0.
%
%   See also MAPTRIMS, MAPTRIML.

% Copyright 1996-2008 The MathWorks, Inc.
% $Revision: 1.1.6.4 $  $Date: 2008/06/16 16:47:52 $

%  Make sure latlim and lonlim are 2-vectors.
if (numel(latlim) ~= 2) || (numel(lonlim) ~= 2)
    eid = sprintf('%s:%s:wrongSizedLimits',getcomp,mfilename); 
    error(eid,'Lat and lon limit inputs must be two-element vectors.')
end

% Condition the longitude limits.
lonlim = conditionLonlim(lonlim);

% Convert inputs to radians and trimPolygonToQuadrangle.  Set inc
% to 2*pi to avoid interpolating additional points along the edges of
% the quadrangle.
inc = 2*pi;
[lat, lon] = trimPolygonToQuadrangle(...
    lat*pi/180, lon*pi/180, latlim*pi/180, lonlim*pi/180, inc);

% Convert back to degrees.
lat = lat*180/pi;
lon = lon*180/pi;
