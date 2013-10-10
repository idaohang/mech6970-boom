function R = maprasterref(varargin)
%maprasterref Construct spatialref.MapRasterReference object
%
%   R = maprasterref() constructs a spatialref.MapRasterReference object with
%   default property values.
%
%   R = maprasterref(Prop1Name, Prop1Value, Prop2Name, Prop2Value, ...)
%   accepts a list of name-value pairs that are used to assign selected
%   properties when initializing a spatialref.MapRasterReference object.
%   You can include any of the following properties, overriding their
%   default values as needed:
%
%       XLimWorld (default value: [0.5 2.5])
%       YLimWorld (default value: [0.5 2.5])
%       RasterSize (default value: [2 2])
%       RasterInterpretation (default value: 'cells')
%       ColumnsStartFrom (default value: 'south')
%       RowsStartFrom (default value: 'west')
%
%   Alternatively, you may omit any or all properties when constructing
%   your spatialref.MapRasterReference object. Then you can customize the
%   result by resetting properties from this list one at a time.
%
%   This name-value syntax always results in an object with a
%   'rectilinear' TransformationType. If your image is rotated with
%   respect to the world coordinate axes, you need an object with a
%   TransformationType of 'affine'. You can obtain such an object
%   directly from the spatialref.MapRasterReference constructor. 
%   Alternately, you can provide an appropriate world file matrix as input,
%   as shown in the following syntax. You cannot do it by resetting
%   properties of an existing rectilinear spatialref.MapRasterReference
%   object.
%
%   R = maprasterref(W, rasterSize, rasterInterpretation) constructs a
%   spatialref.MapRasterReference object with the specified raster size
%   and interpretation properties, and with remaining properties defined
%   by a 2-by-3 world file matrix, W. The rasterInterpretation input is
%   optional, can equal either 'cells' or 'postings', and has a default
%   value of 'cells'.
%
%   Example 1
%   ---------
%   Construct a referencing object for an 1000-by-2000 image with
%   square, 1/2 meter pixels referenced to a planar map coordinate
%   system (the "world" system). The X-limits in the world system are
%   207000 and 208000. The Y-limits are 912500 and 913000. The image
%   follows the popular convention in which world X increases from
%   column to column and world Y decreases from row to row.
%
%   % Override the default MATLAB display format. This is not strictly
%   % required, but tends to produces the most readable displays.
%   format short g
%
%   % Construct a spatialref.MapRasterReference object.
%   R = maprasterref('RasterSize', [1000 2000], ...
%         'YLimWorld', [912500 913000], 'ColumnsStartFrom','north', ...
%         'XLimWorld', [207000 208000])
%
%   Example 2
%   ---------
%   % Repeat Example 1 with a different strategy: Create a default
%   % object and then modify that object's property settings as needed.
%   R = maprasterref;
%   R.XLimWorld = [207000 208000];
%   R.YLimWorld = [912500 913000];
%   R.ColumnsStartFrom = 'north';
%   R.RasterSize = [1000 2000]
%
%   Example 3
%   ---------
%   % Repeat Example 1 again, this time using a world file matrix.
%   W = [0.5   0.0   207000.25; ...
%        0.0  -0.5   912999.75];
%   rasterSize = [1000 2000];
%   R = maprasterref(W, rasterSize)
%
%   See also georasterref, spatialref.MapRasterReference

% Copyright 2010-2011 The MathWorks, Inc.
% $Revision: 1.1.6.2.2.1 $  $Date: 2011/02/05 19:22:27 $

if nargin == 0 || ischar(varargin{1})
    % Construct a default object.
    R = spatialref.MapRasterReference();
    
    % Allow the following properties to be set.
    validPropertyNames = {'RasterSize', 'RasterInterpretation', ...
        'XLimWorld', 'YLimWorld', 'ColumnsStartFrom', 'RowsStartFrom'};
    
    % Set each property found in the varargin list.
    R = setSpatialReferencingProperties(R, ...
        varargin, validPropertyNames, 'maprasterref');
else
    W = varargin{1};  % We know there's at least one input.
    W = validateWorldFileMatrix(W, 'maprasterref', 'worldFileMatrix', 1);
    
    [rasterSize, rasterInterpretation] ...
        = parseRasterSizeAndInterpretation(varargin(2:end), 'maprasterref');
    
    J = W(:,1:2);
    
    if strcmp(rasterInterpretation, 'cells')
        firstCornerX = firstCorner(W(1,3), -(J(1,1) + J(1,2)));
        firstCornerY = firstCorner(W(2,3), -(J(2,1) + J(2,2)));
    else
        firstCornerX = W(1,3);
        firstCornerY = W(2,3);
    end
    
    R = spatialref.MapRasterReference(rasterSize, rasterInterpretation, ...
        '', firstCornerX, firstCornerY, J, [1 1; 1 1]);
end
