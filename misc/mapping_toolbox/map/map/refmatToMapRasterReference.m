function R = refmatToMapRasterReference( ...
    refmat, rasterSize, func_name, var_name, arg_pos)
%refmatToMapRasterReference Referencing matrix to MapRasterReference object
%
%   R = refmatToMapRasterReference(REFMAT, rasterSize) constructs a
%   spatialref.MapRasterReference object, R, from a referencing matrix,
%   REFMAT, and size vector, rasterSize.
%
%   rasterSize is a size vector [M N ...] specifying the number of rows
%   (M) and columns (N) in the raster or image to be associated with the
%   MapRasterReference object, R. For convenience, rasterSize may be a
%   row vector with more than two elements. This flexibility enables you
%   to specify the size in the following way:
%
%       R = refmatToMapRasterReference(refmat, size(RGB))
%
%   where RGB is M-by-N-by-3. However, in such cases, only the first
%   two elements of the size vector are actually used. The higher
%   (non-spatial) dimensions are ignored.
%
%   R = refmatToMapRasterReference(REFMAT, rasterSize, FUNC_NAME, VAR_NAME, ARG_POS)
%   uses up to three optional arguments to provide additional
%   information. This information is used to construct error messages if
%   either the REFMAT or rasterSize inputs turn out to be invalid. Thus,
%   you can use refmatToMapRasterReference for both validating and
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
%      refmatToMapRasterReference includes this information in the
%      formatted error message.
%
%   R = refmatToMapRasterReference(Rin, rasterSize, ...), where Rin is
%   a spatialref.MapRasterReference object, verifies that Rin.RasterSize
%   is consistent with rasterSize, then copies Rin to R.
%
%   Example
%   -------
%   % Import a referencing matrix from a world file for a 2000-by-2000
%   % orthoimage referenced to the Massachusetts State Plane Mainland
%   % coordinate system.
%   refmat = worldfileread('concord_ortho_e.tfw')
%
%   % Import the corresponding TIFF image and use its size to help
%   % convert the referencing matrix to a referencing object.
%   [X, cmap] = imread('concord_ortho_e.tif');
%   R = refmatToMapRasterReference(refmat, size(X))
%
%   % Use the mapbbox function to obtain the map limits independently of
%   % the referencing object.
%   bbox = mapbbox(refmat, size(X))
%   xLimWorld = bbox(:,1)';  % Transpose the first column
%   yLimWorld = bbox(:,2)';  % Transpose the second column
%
%   % Construct a referencing object directly, for comparison.
%   maprasterref('RasterSize', size(X), 'ColumnsStartFrom', 'north', ...
%       'XLimWorld', xLimWorld, 'YLimWorld', yLimWorld)
%
%   See also maprasterref, refmatToGeoRasterReference

% Copyright 2010-2011 The MathWorks, Inc.
% $Revision: 1.1.6.2.2.1 $  $Date: 2011/02/05 19:22:31 $

if nargin < 3
    func_name = 'refmatToMapRasterReference';
end

if nargin < 4
    var_name = 'REFMAT';
end

if nargin < 5
    arg_pos = 1;
end

% Validate first input (refmat). It must be a 3-by-2 matrix of real-valued
% finite doubles, or a spatialref.MapRasterReference object.
internal.map.validateRasterReference(refmat, ...
    {'spatialref.MapRasterReference'}, func_name, var_name, arg_pos)

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
        'A %s object was provided to function %s, but its %s value is not consistent with the raster size vector.', ...
        class(R), func_name, 'RasterSize')
else
    % Assume 'cells'.
    rasterInterpretation = 'cells';
    
    isRectilinear = (refmat(1,1) == 0 && refmat(2,2) == 0);
    if isRectilinear
        % Obtain defining parameters from referencing matrix.
        deltaX = refmat(2,1);
        deltaY = refmat(1,2);
        
        firstCornerX = firstCorner(refmat(3,1), deltaX);
        firstCornerY = firstCorner(refmat(3,2), deltaY);
        
        J = [deltaX 0; 0 deltaY];
    else
        % Permute the first two rows of refmat to obtain the Jacobian matrix.
        J = refmat([2 1; 5 4]);
        
        % The bottom row of refmat establishes a reference point at (0, 0)
        % in intrinsic coordinates, so take a step of 0.5 in both X and Y to
        % reach the outer corner at (0.5, 0.5).
        firstCornerX = refmat(3,1) + (J(1,1) + J(1,2)) / 2;
        firstCornerY = refmat(3,2) + (J(2,1) + J(2,2)) / 2;
    end
    
    % Construct referencing object.
    R = spatialref.MapRasterReference( rasterSize, rasterInterpretation, ...
        '', firstCornerX, firstCornerY, J, [1 1; 1 1]);
end
