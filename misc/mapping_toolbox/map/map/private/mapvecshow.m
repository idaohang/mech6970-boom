function varargout = mapvecshow(x, y, varargin)
%MAPVECSHOW Display x-y map vectors with given geometry
%
%   MAPVECSHOW(X,Y) or 
%   MAPVECSHOW(X,Y, ..., 'DisplayType', DISPLAYTYPE, ...) displays the
%   equal length coordinate vectors, X and Y.  X and Y may contain embedded
%   NaNs, delimiting coordinates of points, lines or polygons. DISPLAYTYPE
%   can be 'point', 'multipoint', 'line', or 'polygon' and defaults to
%   'line'.
%
%   MAPVECSHOW(..., PARAM1, VAL1, PARAM2, VAL2, ...) specifies
%   parameter/value pairs that set MATLAB graphics properties. Parameter
%   names can be abbreviated and are case-insensitive.
%
%   H = MAPVECSHOW(...) returns a handle to a MATLAB graphics object.
%
%   Example 
%   -------
%   % Display the Boston roads as black lines.
%   roads = shaperead('boston_roads.shp');
%   figure
%   mapvecshow([roads.X], [roads.Y], 'Color', 'black');
%
%   See also GEOVECSHOW, MAPSHOW, MAPSTRUCTSHOW.

% Copyright 2006-2009 The MathWorks, Inc.
% $Revision: 1.1.6.2 $  $Date: 2009/11/09 16:26:08 $

% Verify the coordinate arrays.
checkxy(x, y, 'MAPSHOW', 'X', 'Y', 1, 2);

% Find the geometry.
% Delete 'DisplayType' parameter/value pair from varargin (if present).
% If DisplayType is not set, then draw a line.
default = 'line';
[geometry, varargin] = ...
   internal.map.findNameValuePair('DisplayType', default, varargin{:});

% Split HG property value pairs into two groups, depending on whether or
% not a (case-insensitive) prefix of 'Default' is included in the
% property name.
[defaultProps, otherProps] = separateDefaults(varargin);
   
% Determine the display function based on the geometry.
fcn = mapvecfcn(geometry, 'mapshow');

% Display the x, y vectors.
h = fcn(x, y, defaultProps{:}, otherProps{:});

% Allow usage without ; to suppress output.
if nargout > 0
   varargout{1} = h;
end
