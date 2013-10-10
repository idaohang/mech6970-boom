function val = ltln2val(Z, R, lat, lon, method)
%LTLN2VAL  Extract data grid values for specified locations
%
%   VAL = LTLN2VAL(Z, R, LAT, LON) interpolates a regular data grid Z at
%   the points specified by vectors of latitude and longitude, LAT and
%   LON.  R can be a spatialref.GeoRasterReference object, a
%   referencing vector, or a referencing matrix.
%
%   If R is a spatialref.GeoRasterReference object, its RasterSize
%   property must be consistent with size(Z).
%
%   If R is a referencing vector, it must be a 1-by-3 with elements:
%
%     [cells/degree northern_latitude_limit western_longitude_limit]
%
%   If R is a referencing matrix, it must be 3-by-2 and transform raster
%   row and column indices to/from geographic coordinates according to:
% 
%                     [lon lat] = [row col 1] * R.
%
%   If R is a referencing matrix, it must define a (non-rotational,
%   non-skewed) relationship in which each column of the data grid falls
%   along a meridian and each row falls along a parallel.
%   Nearest-neighbor interpolation is used by default.  NaN is returned
%   for points outside the grid limits or for which LAT or LON contain
%   NaN.  All angles are in units of degrees.
%
%   VAL = LTLN2VAL(Z, R, LAT, LON, METHOD) accepts a METHOD string
%   to specify the type of interpolation: 'bilinear' for linear
%   interpolation, 'bicubic' for cubic interpolation, or 'nearest' for
%   nearest neighbor interpolation.
%
%   See also FINDM, IMBEDM.

% Copyright 1996-2010 The MathWorks, Inc.
% $Revision: 1.11.4.11 $  $Date: 2010/11/17 11:24:21 $

error(nargchk(4, 5, nargin, 'struct'))
if nargin == 4
    method = 'nearest';
else
    method = validatestring(method, ...
        {'nearest', 'linear', 'cubic', 'bilinear', 'bicubic'}, ...
        'ltln2val', 'METHOD', 5);
    if strcmp(method(1:2), 'bi')
        % Convert 'binear' to 'linear' and 'bicubic' to 'cubic'.
        method(1:2) = [];
    end
end

lat = ignoreComplex(lat, mfilename, 'LAT');
lon = ignoreComplex(lon, mfilename, 'LON');
checklatlon(lat, lon, mfilename, 'LAT', 'LON', 3, 4)
validateattributes(Z,{'numeric','logical'},{'2d','real'},mfilename,'Z',1)

%  If R is already spatial referencing object, validate it. Otherwise
%  convert the input referencing vector or matrix.
R = internal.map.convertToGeoRasterRef( ...
    R, size(Z), 'degrees', 'LTLN2VAL', 'R', 2);

%  Remove NaNs from lat/lon arrays, but keep track of where they were.
nanLatOrLon = isnan(lat) | isnan(lon);
lat(nanLatOrLon) = [];
lon(nanLatOrLon) = [];

% Initialize output to an array of NaN matching the original size of lat
% and lon, then fill in values only for elements corresponding to
% non-NaN lat and lon.
val = NaN(size(nanLatOrLon));
if method(1) == 'n' && strcmp(R.RasterInterpretation,'cells')
    % If the interpolation method is 'nearest', use SETPOSTN.
    val(~nanLatOrLon) = interpNearest(Z, R, lat, lon);
else
    % Use INTERP2.
    val(~nanLatOrLon) = interpGeoRaster(Z, R, lat, lon, method);
end

%-----------------------------------------------------------------------

function val = interpNearest(Z, R, lat, lon)

% Nearest-neighbor interpolation; return a value of NaN for points
% outside the grid limits.  VAL will be the same size as LAT and LON.

% Turn off warning for points outside map limits during call to SETPOSTN
w = warning('off','map:setpostn:pointOutsideLimits');
c = onCleanup(@() warning(w));
[row, col, indxPointOutsideGrid] = setpostn(Z, R, lat, lon);

% Construct logical array matching LAT and LON in size to indicate which
% points are inside or outside the grid limits.
insideLimits = true(size(lat));
insideLimits(indxPointOutsideGrid) = false;

% Output array matches LAT and LON in size and contains values copied
% from Z for points inside the map limits and NaN elsewhere.
indx = row + size(Z,1)*(col - 1);
val = NaN(size(lat));
val(insideLimits) = Z(indx);
