function R = georasterref(varargin)
%georasterref Construct spatialref.GeoRasterReference object
%
%   R = georasterref() constructs a spatialref.GeoRasterReference object with
%   default property values.
%
%   R = georasterref(Prop1Name, Prop1Value, Prop2Name, Prop2Value, ...)
%   accepts a list of name-value pairs that are used to assign selected
%   properties when initializing a spatialref.GeoRasterReference object.
%   You can include any of the following properties, overriding their
%   default values as needed:
%
%       Latlim (default value: [0.5 2.5])
%       Lonlim (default value: [0.5 2.5])
%       RasterSize (default value: [2 2])
%       RasterInterpretation (default value: 'cells')
%       ColumnsStartFrom (default value: 'south')
%       RowsStartFrom (default value: 'west')
%
%   Alternatively, you may omit any or all properties when constructing
%   your spatialref.GeoRasterReference object. Then you can customize the
%   result by resetting properties from this list one at a time.
%
%   R = georasterref(W, rasterSize, rasterInterpretation) constructs a
%   spatialref.GeoRasterReference object with the specified raster size
%   and interpretation properties, and with remaining properties defined
%   by a 2-by-3 world file matrix, W. The rasterInterpretation input is
%   optional, can equal either 'cells' or 'postings', and has a default
%   value of 'cells'.
%
%   Example 1
%   ---------
%   % Construct a referencing object for a global raster
%   % comprising 180-by-360 one-degree cells, with rows that
%   % start at longitude -180, and with the first cell
%   % located in the northwest corner.
%
%   % Override the default MATLAB display format. This is not strictly
%   % required, but tends to produces the most readable displays.
%   format short g
%
%   % Construct a spatialref.GeoRasterReference object.
%   R = georasterref('RasterSize', [180 360], ...
%         'RasterInterpretation', 'cells', ...
%         'Latlim', [-90 90], 'Lonlim', [-180 180], ...
%         'ColumnsStartFrom', 'north')
%
%   Example 2
%   ---------
%   % Construct a referencing object for the DTED Level 0 file that
%   % includes Sagarmatha (Mount Everest). The DTED columns run
%   % from south to north and the first column runs along the
%   % western edge of the (one-degree-by-one-degree) quadrangle,
%   % consistent with the default values for 'ColumnsStartFrom' and
%   % 'RowsStartFrom'.
%   R = georasterref('Latlim', [27 28], 'Lonlim', [86 87], ...
%        'RasterSize', [121 121], 'RasterInterpretation', 'postings')
%
%   Example 3
%   ---------
%   % Repeat Example 2 with a different strategy: Create a default
%   % object and then modify its properties as needed.
%   R = georasterref;
%   R.RasterSize = [121 121];
%   R.RasterInterpretation = 'postings';
%   R.Latlim = [27 28];
%   R.Lonlim = [86 87]
%
%   Example 4
%   ---------
%   % Repeat Example 1 using a world file matrix as input.
%   W = [1    0   -179.5; ...
%        0   -1     89.5];
%   rasterSize = [180 360];
%   rasterInterpretation = 'cells';
%   R = georasterref(W, rasterSize, rasterInterpretation);
%
%   See also maprasterref, spatialref.GeoRasterReference

% Copyright 2010-2011 The MathWorks, Inc.
% $Revision: 1.1.6.3.2.1 $  $Date: 2011/02/05 19:22:25 $

if nargin == 0 || ischar(varargin{1})
    % Construct a default object.
    R = spatialref.GeoRasterReference();
    
    % Allow the following properties to be set.
    validPropertyNames = {'RasterSize', 'RasterInterpretation', ...
        'Latlim', 'Lonlim', 'ColumnsStartFrom', 'RowsStartFrom'};
    
    % Set each property found in the varargin list.
    R = setSpatialReferencingProperties(R, ...
        varargin, validPropertyNames, 'georasterref');
else
    W = varargin{1};  % We know there's at least one input.
    W = validateWorldFileMatrix(W, 'georasterref', 'worldFileMatrix', 1);
    
    assert(W(1,2) == 0 && W(2,1) == 0, ...
        'map:validate:expectedRectilinearWorldFileMatrix', ...
        'Function %s expected input number %d, %s, to define a rectilinear transformation between intrinsic and geographic coordinates.', ...
        'georasterref', 1, 'worldFileMatrix')
    
    [rasterSize, rasterInterpretation] ...
        = parseRasterSizeAndInterpretation(varargin(2:end), 'georasterref');
    
    deltaLon = W(1,1);
    deltaLat = W(2,2);
    
    if strcmp(rasterInterpretation, 'cells')
        firstCornerLon = firstCorner(W(1,3), -deltaLon);
        firstCornerLat = firstCorner(W(2,3), -deltaLat);
    else
        firstCornerLon = W(1,3);
        firstCornerLat = W(2,3);
    end
    
    R = spatialref.GeoRasterReference(rasterSize, rasterInterpretation, ...
        'degrees', firstCornerLat, firstCornerLon, deltaLat, 1, deltaLon, 1);
end
