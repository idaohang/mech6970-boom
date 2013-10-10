function [lat, lon] = projinv(proj, x, y)
%PROJINV Inverse map projection using PROJ.4 library
%
%   [LAT, LON] = PROJINV(PROJ, X, Y) returns the latitude and longitude
%   values from the inverse projection.  PROJ is a structure defining the
%   map projection.  PROJ may be a map projection MSTRUCT or a GeoTIFF INFO
%   structure.  X and Y are map coordinate arrays.  For a complete list of
%   GeoTIFF info and map projection structures that may be used with
%   PROJINV, see PROJLIST.
%
%   Example - Overlay Boston roads on top of 'boston_ovr.jpg' 
%             in a Mercator projection.
%   -------
%   % Import the Boston roads from the shapefile.
%   roads = shaperead('boston_roads.shp');
%
%   % Obtain the projection structure from 'boston.tif'.
%   proj = geotiffinfo('boston.tif');
%
%   % As shown by the UOMLength field of the projection structure, the
%   % units of length in the projected coordinate system is 
%   % 'US Survey Feet'. Coordinates in the roads shapefile are in meters.
%   % The road coordinates must be converted to the projection's length 
%   % unit.
%   x = [roads.X] * unitsratio('sf','meter');
%   y = [roads.Y] * unitsratio('sf','meter');
%  
%   % Convert the coordinates of the roads to latitude and longitude.
%   [roadsLat, roadsLon] = projinv(proj, x, y);
%
%   % Read the boston_ovr.jpg image and worldfile.
%   % Includes material (c) GeoEye, all rights reserved.
%   RGB = imread('boston_ovr.jpg');
%   R = worldfileread(getworldfilename('boston_ovr.jpg'));
%
%   % Obtain stateline boundary of Massachusetts.
%   S = shaperead('usastatehi', 'UseGeoCoords', true, ...
%       'Selector',{@(name) strcmpi(name,'Massachusetts'), 'Name'});
%
%   % Open a figure with a Mercator projection.
%   figure
%   axesm('mercator')
%
%   % Display the stateline boundary, image, and roads.
%   geoshow(S.Lat, S.Lon, 'Color','red')
%   geoshow(RGB, R)
%   geoshow(roadsLat, roadsLon, 'Color', 'green')
%
%   % Set the map boundary to the image's northern, western, and southern
%   % limits, and the eastern limit of the stateline within the image
%   % latitude boundaries.
%   [lon, lat] = mapoutline(R, size(RGB(:,:,1)));
%   ltvals = find((S.Lat>=min(lat(:))) & (S.Lat<=max(lat(:))));
%   setm(gca,'maplonlimit',[min(lon(:)) max(S.Lon(ltvals))], ...
%            'maplatlimit',[min(lat(:)) max(lat(:))])
%   tightmap
%   
%   See also GEOTIFFINFO, MFWDTRAN, MINVTRAN, PROJFWD, PROJLIST.

% Copyright 1996-2008 The MathWorks, Inc.
% $Revision: 1.1.6.9 $  $Date: 2008/05/14 22:01:48 $

% Check the input arguments
checknargin(3,3,nargin,mfilename);
checkinput(x, {'numeric'}, {'real'}, mfilename, 'X', 2);
checkinput(y, {'numeric'}, {'real'}, mfilename, 'Y', 3);
if numel(x) ~= numel(y)
   eid = sprintf('%s:%s:invalidLength', getcomp, mfilename);
   msg = sprintf('The number of elements in X and Y must be equal.');
   error(eid, '%s',msg);
end

% Inverse transform the X and Y points.
[lat, lon] = projaccess('inv', proj, x, y);

