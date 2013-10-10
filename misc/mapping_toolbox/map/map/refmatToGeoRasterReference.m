function R = refmatToGeoRasterReference( ...
    refmat, rasterSize, func_name, var_name, arg_pos)
%refmatToGeoRasterReference Referencing matrix to GeoRasterReference object
%
%   R = refmatToGeoRasterReference(REFMAT, rasterSize) constructs a
%   spatialref.GeoRasterReference object, R, from a referencing matrix,
%   REFMAT, and size vector, rasterSize. REFMAT may be any valid
%   referencing matrix subject to the two following constraints. First, the
%   matrix must lead to valid latitude and longitude limits when combined
%   with rasterSize. Second, the matrix columns and rows must be aligned
%   with meridians and parallels, respectively.
%
%   rasterSize is a size vector [M N ...] specifying the number of rows
%   (M) and columns (N) in the raster or image to be associated with the
%   GeoRasterReference object, R. For convenience, rasterSize may be a
%   row vector with more than two elements. This flexibility enables you
%   to specify the size in the following way:
%
%       R = refmatToGeoRasterReference(refmat, size(RGB))
%
%   where RGB is M-by-N-by-3. However, in such cases, only the first
%   two elements of the size vector are actually used. The higher
%   (non-spatial) dimensions are ignored.
%
%   R = refmatToGeoRasterReference(REFMAT, rasterSize, FUNC_NAME, VAR_NAME, ARG_POS)
%   uses up to three optional arguments to provide additional
%   information. This information is used to construct error messages if
%   either the REFMAT or rasterSize inputs turn out to be invalid. Thus,
%   you can use refmatToGeoRasterReference for both validating and
%   converting a referencing matrix. The optional inputs work just like
%   their counterparts in the MATLAB function VALIDATEATTRIBUTES:
%
%      FUNC_NAME is a string that specifies the name used in the
%      formatted error message to identify the function checking the
%      input.
%
%      VAR_NAME is a string that specifies the name used in the
%      formatted error message to identify the referencing matrix.
%
%      ARG_POS is a positive integer that indicates the position of the
%      referencing matrix checked in the function argument list.
%      refmatToGeoRasterReference includes this information in the
%      formatted error message.
%
%   R = refmatToGeoRasterReference(Rin, rasterSize, ...), where Rin is
%   a spatialref.GeoRasterReference object, verifies that Rin.RasterSize
%   is consistent with rasterSize, then copies Rin to R.
%
%   Example
%   -------
%   % Construct a referencing matrix for a regular grid that covers the
%   % entire globe with 1-degree cells.
%   rasterSize = [180 360];
%   refmat = makerefmat( ...
%       'RasterSize', rasterSize, 'Latlim', [-90 90], 'Lonlim', [0 360])
%
%   % Convert to a spatialref.GeoRasterReference object
%   R = refmatToGeoRasterReference(refmat, rasterSize)
%
%   % For comparison, construct a referencing object directly.
%   georasterref( ...
%       'RasterSize', rasterSize, 'Latlim', [-90 90], 'Lonlim', [0 360])
%
%   See also georasterref, refvecToGeoRasterReference

% Copyright 2010-2011 The MathWorks, Inc.
% $Revision: 1.1.6.1.2.1 $  $Date: 2011/02/05 19:22:29 $

if nargin < 3
    func_name = 'refmatToGeoRasterReference';
end

if nargin < 4
    var_name = 'REFMAT';
end

if nargin < 5
    arg_pos = 1;
end

% Validate first input (refmat). It must be a 3-by-2 matrix of real-valued
% finite doubles, or a spatialref.GeoRasterReference object.
internal.map.validateRasterReference(refmat, ...
    {'spatialref.GeoRasterReference'}, func_name, var_name, arg_pos)

%Validate raster size.
validateattributes(rasterSize, ...
    {'double'}, {'row','positive','integer'}, func_name)

assert(numel(rasterSize) >= 2, ...
    'map:convertspatialref:invalidRasterSize', ...
    'The raster size vector provided to function %s must have at least two elements.', ...
    func_name)

if isobject(refmat)
    R = refmat;
    
    assert(isequal(R.RasterSize, rasterSize(1,1:2)), ...
        'map:convertspatialref:inconsistentRasterSize', ...
        'A %s object was provided to function %s., but its % value is not consistent with the raster size vector.', ...
        class(R), func_name, 'RasterSize')
else
    assert(refmat(1,1) == 0 && refmat(2,2) == 0, ...
        'map:convertspatialref:skewOrRotation', ...
        ['The referencing matrix supplied to function %s specifies', ...
        ' that the associated raster is rotated or skewed with', ...
        ' respect to the latitude/longitude system.  Function %s', ...
        ' does not support this geometry.'], func_name, func_name)
    
    % Assume 'cells' and 'degrees'.
    rasterInterpretation = 'cells';
    angleUnits = 'degrees';
    
    % Compute other defining parameters from referencing matrix.
    deltaLat = refmat(1,2);
    deltaLon = refmat(2,1);
    
    firstCornerLat = firstCorner(refmat(3,2), deltaLat);
    firstCornerLon = firstCorner(refmat(3,1), deltaLon);
    
    % Invoke constructor and interpret errors in the case of inconsistent
    % combinations of property values.
    R = constructGeoRasterReference( ...
        rasterSize, rasterInterpretation, angleUnits, ...
        firstCornerLat, firstCornerLon, deltaLat, 1, deltaLon, 1);
end
