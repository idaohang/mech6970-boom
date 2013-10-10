function [lat, lon] = pix2latlon(R, row, col)
%PIX2LATLON Convert pixel coordinates to latitude-longitude coordinates
%
%   [LAT, LON] = PIX2LATLON(R,ROW,COL) calculates latitude-longitude
%   coordinates LAT, LON from pixel coordinates ROW, COL.  R is either a
%   3-by-2 referencing matrix that transforms intrinsic pixel coordinates
%   to geographic coordinates, or a spatialref.GeoRasterReference object.
%   ROW and COL are vectors or arrays of matching size.  The outputs LAT
%   and LON have the same size as ROW and COL.
%
%   Example 
%   -------
%      % Find the latitude and longitude of the upper left and lower right 
%      % outer corners of a 2-by-2 degree gridded data set.
%      R = makerefmat(1, 89, 2, 2);
%      [UL_lat, UL_lon] = pix2latlon(R, .5, .5)
%      [LR_lat, LR_lon] = pix2latlon(R, 90.5, 180.5)
%
%   See also LATLON2PIX, MAKEREFMAT, PIX2MAP.

% Copyright 1996-2010 The MathWorks, Inc.
% $Revision: 1.1.8.6 $  $Date: 2010/10/11 14:47:23 $ 

error(nargchk(3,3,nargin,'struct'))

% Validate referencing matrix or GeoRasterReference object.
internal.map.validateRasterReference(R, ...
    {'spatialref.GeoRasterReference'}, 'pix2latlon', 'R', 1)

if isobject(R)
    [lat, lon] = R.intrinsicToGeographic(col, row);
else
    [lon, lat] = pix2map(R, row, col);
end
