function varargout = geovecshow(lat, lon, varargin)
%GEOVECSHOW Display geographic vectors with given geometry
%
%   GEOVECSHOW(LAT, LON) or
%   GEOVECSHOW(LAT, LON, ..., 'DisplayType', DISPLAYTYPE, ...) displays the
%   equal length coordinate vectors, LAT and LON.  LAT and LON may contain
%   embedded NaNs, delimiting coordinates of points, lines or polygons.
%   DISPLAYTYPE can be 'point', 'multipoint', 'line', or 'polygon' and
%   defaults to 'line'.
%
%   GEOVECSHOW(..., PARAM1, VAL1, PARAM2, VAL2, ...) specifies
%   parameter/value pairs that set MATLAB graphics properties. Parameter
%   names can be abbreviated and are case-insensitive.
%
%   H = GEOVECSHOW(...) returns a handle to a MATLAB graphics object.
%
%   See also GEOSHOW, GEOSTRUCTSHOW.
%
%   Example 1 
%   ---------
%   % Display world coastlines, using a Plate Carree projection.
%   load coast
%   figure
%   geovecshow(lat, long);
%
%   Example 2 
%   ---------
%   % Create a worldmap of North America and draw land areas as black lines
%   land = shaperead('landareas.shp','UseGeoCoords', true);
%   figure
%   worldmap('na');
%   geovecshow([land.Lat], [land.Lon], 'Color', 'black')

% Copyright 2006-2009 The MathWorks, Inc.
% $Revision: 1.1.6.3 $  $Date: 2009/11/09 16:26:03 $

% Verify the coordinate arrays.
checklatlon(lat, lon, 'GEOSHOW', 'LAT', 'LON', 1, 2);

% Find the geometry.
% Delete 'DisplayType' parameter/value pair from varargin (if present).
% If DisplayType is not set, then draw a line.
default = 'line';
[geometry, varargin] = ...
   internal.map.findNameValuePair('DisplayType', default, varargin{:});

% Obtain the projection structure.
mstruct = getProjection(varargin{:});

% Split HG property value pairs into two groups, depending on whether or
% not a (case-insensitive) prefix of 'Default' is included in the property
% name.
[defaultProps, otherProps] = separateDefaults(varargin);

% Switch display and projection operations based on Geometry.
objectType = ['geo' lower(geometry)];

if strcmpi(mstruct.mapprojection,'globe')
    % Treat globe axes separately because of the third dimension.
    height = [];
    h = globevec(mstruct, lat, lon, height, ...
        objectType, defaultProps{:}, otherProps{:});
else
    % Project and display in 2-D map coordinates.
    fcn = mapvecfcn(geometry, 'geoshow');
    h = geovec(mstruct, lat, lon,  ...
        objectType, fcn, defaultProps{:}, otherProps{:});
end

% Allow usage without ; to suppress output.
if nargout > 0
   varargout{1} = h;
end
