function [x, y] = projfwd(proj, lat, lon)
%PROJFWD Forward map projection using PROJ.4 library
%
%   [X, Y] = PROJFWD(PROJ, LAT, LON) returns the X and Y map coordinates
%   from the forward projection.  PROJ is a structure defining the map
%   projection.  PROJ may be a map projection MSTRUCT or a GeoTIFF INFO
%   structure.  LAT and LON are arrays of latitude and longitude
%   coordinates. For a complete list of GeoTIFF info and map projection
%   structures that may be used with PROJFWD, see PROJLIST.
%
%   Example 1
%   ---------
%   % Overlay landarea boundary on 'boston.tif'.
%   % Includes material (c) GeoEye, all rights reserved.
%
%   % Obtain stateline boundary of Massachusetts.
%   S = shaperead('usastatehi', 'UseGeoCoords', true, ...
%       'Selector',{@(name) strcmpi(name,'Massachusetts'), 'Name'});
%
%   % Obtain the projection structure.
%   proj = geotiffinfo('boston.tif');
%
%   % Project the stateline boundary.
%   lat = [S.Lat];
%   lon = [S.Lon];
%   [x, y] = projfwd(proj, lat, lon);
%
%   % Read the 'boston.tif' image.
%   [RGB, R] = geotiffread('boston.tif');
%
%   % Display the image.
%   figure
%   mapshow(RGB, R)
%   xlabel('MA Mainland State Plane easting, survey feet')
%   ylabel('MA Mainland State Plane northing, survey feet')
%
%   % Overlay the stateline boundary.
%   hold on
%   mapshow(gca, x, y,'Color','black','LineWidth',2.0)
%
%   % Set the map boundary to show a little more detail.
%   set(gca,'XLim', [ 645000,  895000], ...
%           'YLIm', [2865000, 3040000]);
%
%   See also GEOTIFFINFO, MFWDTRAN, MINVTRAN, PROJINV, PROJLIST.

% Copyright 1996-2010 The MathWorks, Inc.
% $Revision: 1.1.6.11 $  $Date: 2010/11/17 11:25:20 $

% Check the input arguments
checknargin(3,3,nargin,mfilename);
checkinput(lat, {'numeric'}, {'real'}, mfilename, 'LAT', 2);
checkinput(lon, {'numeric'}, {'real'}, mfilename, 'LON', 3);
if numel(lat) ~= numel(lon)
   eid = sprintf('%s:%s:invalidLength', getcomp, mfilename);
   msg = sprintf('The number of elements in LAT and LON must be equal.');
   error(eid, '%s',msg);
end

% Project the latitude and longitude points.
[x,y] = projaccess('fwd', proj, lat, lon);
