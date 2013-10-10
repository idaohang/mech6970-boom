function [outputR, rows, cols] =  selectGeoRasterSubsamplingLimits( ...
    inputR, latlim, lonlim, rowSampleFactor, columnSampleFactor, ...
    columnsRunSouthToNorth, rowsRunWestToEast)
%selectGeoRasterSubsamplingLimits Select indices for cell-oriented raster
%subsampling.
%
%  [outputR, rows, cols] = selectGeoRasterSubsamplingLimits(inputR, ...
%  latlim, lonlim, rowSampleFactor, columnSampleFactor, ...
%  columnsRunSouthToNorth, rowsRunWestToEast) selects the row and column
%  indices to use when subsampling a cell-oriented raster grid, given
%  user-specified limits and sample factor.
%
%  Inputs:
%     inputR - a GeoRasterReference object defined for the input raster grid
%              with RasterInterpretation equal to 'cells'
%           
%     latlim - a two-element latitude limits vector of the form 
%              [southern_limit northern_limit]
%
%     lonlim - a two-element longitude limits vector of the form 
%              [western_limit eastern_limit]
%
%     rowSampleFactor - scalar subsampling factor to be applied in the 
%              row/latitude dimension; must be a positive integer 
%
%     columnSamplefactor - scalar subsampling factor to be applied in the 
%              column/longitude dimension; must be a positive integer 
%
%     columnsRunSouthToNorth - scalar logical; true if and only if columns 
%              are to run from south to north in the output grid 
%
%     rowsRunWestToEast - scalar logical; true if and only if rows are to 
%              run from west to east in the output grid 
%
%  Outputs:
%     outputR - a GeoRasterReference object defined for the output raster grid 
%               with RasterInterpretation equal to 'cells'
%
%     rows    - a row vector indicating which input raster rows are to be 
%               sampled to construct the output raster 
%
%     cols    - a row vector indicating which input raster columns are to
%               be sampled to construct the output raster 
%
%   See also spatialref.GeoRasterReference

% Copyright 2009-2010 The MathWorks, Inc.
% $Revision: 1.1.6.4 $  $Date: 2010/09/24 14:33:37 $

% Validate the requested quadrangle.
[latlim, lonlimInt] = intersectgeoquad( ...
    inputR.Latlim, inputR.Lonlim, latlim, lonlim);
if isempty(latlim)
    error('map:selectGeoRasterSubsamplingLimits:noIntersection', ...
        'The requested limits fail to intersect the input raster grid limits.');
elseif numel(lonlimInt) == 4
    error('map:selectGeoRasterSubsamplingLimits:twoIntersections', ...
        ['Expected the requested longitude interval, %s, to intersect ', ...
        'the limits of the input raster grid exactly once.'], 'LONLIM');
end
if ~isequal(lonlimInt, [-180 180])
    lonlim = lonlimInt;
end

% Validate sample factors. If the sample factor is greater than the size of
% the input raster, then adjust the sample factor.
validateattributes(rowSampleFactor, {'numeric'}, ...
    {'scalar', 'integer', 'positive', 'finite'}, ...
    mfilename, 'rowSampleFactor');
validateattributes(columnSampleFactor, {'numeric'}, ...
    {'scalar', 'integer', 'positive', 'finite'}, ...
    mfilename, 'columnSampleFactor');
rowSampleFactor = min([rowSampleFactor, inputR.RasterSize(1)]);
columnSampleFactor = min([columnSampleFactor, inputR.RasterSize(2)]);

% Process each of the dimensions independently and calculate the output
% row/column indices and latitude/longitude limit.

% Latitude
% Compute the intrinsic coordinates of the latitude limits. If the limits
% run north-to-south, flip them, since the selectSubsamplingLimits function
% require ascending values. 
intrinsicY = inputR.latitudeToIntrinsicY(latlim);
if diff(intrinsicY) < 0
    intrinsicY = intrinsicY([2 1]);
end

% The number of rows in the input raster.
numRows = inputR.RasterSize(1);

% Calculate the intrinsic beginning and ending coordinates and the first
% and last row indices.
[a, b, firstRow, lastRow] = selectSubsamplingLimits( ...
    numRows, rowSampleFactor, intrinsicY(1), intrinsicY(2));

% Calculate the firstCornerLat based on the returned intrinsic coordinates
% [a,b]. Flip the rows if the user requested direction does not match the
% input data grid direction.
inputColumnsRunSouthToNorth = (inputR.DeltaLatNumerator > 0);
sameDirection = (inputColumnsRunSouthToNorth == columnsRunSouthToNorth);
if sameDirection
    rows = firstRow:rowSampleFactor:lastRow;
    if isempty(rows)
        rows = firstRow;
    end
    firstCornerLat = inputR.intrinsicYToLatitude(a);
    deltaLatNumerator = rowSampleFactor * inputR.DeltaLatNumerator;
else
    rows = lastRow:-rowSampleFactor:firstRow;
    if isempty(rows)
        rows = lastRow;
    end
    firstCornerLat = inputR.intrinsicYToLatitude(b);
    deltaLatNumerator = -rowSampleFactor * inputR.DeltaLatNumerator;
end
deltaLatDenominator = inputR.DeltaLatDenominator;

% Longitude
% Compute the intrinsic coordinates of the longitude limits. If the limits
% run east-to-west, flip them, since the selectSubsamplingLimits function
% require ascending values.
intrinsicX = inputR.longitudeToIntrinsicX(lonlim);
if diff(intrinsicX) < 0
    intrinsicX = intrinsicX([2 1]);
end

% The number of columns in the input raster.
numCols = inputR.RasterSize(2);

% Calculate the intrinsic beginning and ending coordinates and the first
% and last row indices.
[a, b, firstCol, lastCol] = selectSubsamplingLimits( ...
    numCols, columnSampleFactor, intrinsicX(1), intrinsicX(2));

% Calculate the firstCornerLon based on the returned intrinsic coordinates
% [a,b]. Flip the columns if the user requested direction does not match
% the input data grid direction.
inputRowsRunWestToEast = (inputR.DeltaLonNumerator > 0);
sameDirection = (inputRowsRunWestToEast == rowsRunWestToEast);
if sameDirection
    cols = firstCol:columnSampleFactor:lastCol;
    if isempty(cols)
        cols = firstCol;
    end
    firstCornerLon = inputR.intrinsicXToLongitude(a);
    deltaLonNumerator = columnSampleFactor * inputR.DeltaLonNumerator;
else
    cols = lastCol:-columnSampleFactor:firstCol;
    if isempty(cols)
        cols = lastCol;
    end
    firstCornerLon = inputR.intrinsicXToLongitude(b);
    deltaLonNumerator = -columnSampleFactor * inputR.DeltaLonNumerator;
end
deltaLonDenominator = inputR.DeltaLonDenominator;

% Create a GeoRasterReference object for the output data grid.
rasterSize = [numel(rows) numel(cols)];

rasterInterpretation = inputR.RasterInterpretation;
angleUnits = inputR.AngleUnits;

outputR = spatialref.GeoRasterReference( ...
    rasterSize, rasterInterpretation, angleUnits, ...
    firstCornerLat, firstCornerLon, ...
    deltaLatNumerator, deltaLatDenominator, ...
    deltaLonNumerator, deltaLonDenominator);

%--------------------------------------------------------------------------

function [a, b, first, last] = selectSubsamplingLimits(n, f, a, b)
% Select the first and last row or column indices to use when
% subsampling a cell-oriented raster grid, given user-specified limits,
% and revise the limits to match a set of output cells with cell extent
% exactly equal to f * (input cell extent).
% 
%   n - number of rows or columns in the input grid
%   f - subsampling factor (positive integer; 1 <= f <= n)
%   a - lower limit of requested data interval, in intrinsic coordinates
%   b - upper limit of requested data interval, in intrinsic coordinates
%
% Intrinsic coordinates are defined such that center of the k-th cell
% falls at the value x = k on the real line. The limits of the input
% raster are thus 0.5 and 0.5 + n.
%
% a and b are real numbers such that 0.5 <= a, a < b, and b < 0.5 + n.
% On output, a and b are adjusted to snap to the input grid while
% respecting other constraints as well.

validateattributes(a, {'double'}, {'scalar','real','finite','>=',0.5}, ...
    mfilename, 'a')
validateattributes(b, {'double'}, {'scalar','real','finite'}, ...
    mfilename, 'b')
validateattributes(n - b, {'double'}, {'>=',-0.5},  ...
    mfilename, 'n - b')
validateattributes(b - a, {'double'}, {'nonnegative'}, ...
  mfilename,'b - a')

% Step 1 -- Determine the number of output cells

% Candidate for number of output cells
m = ceil((b - a) / f);

while f*m > n
    % We don't have enough input cells available to create m output
    % cells, so decrease the number of output cells.
    m = m - 1;
end

% Step 2 -- Estimate the indices of the first and last input cells to be
% covered by the output cells, by starting from an integer close to the
% average of a and b, and going (approximately, in the case of even f*m)
% an equal distance to either side. Note that this definition of first
% and last differs from the specification given in the help above.
% That's OK, we'll use these intermediate values, then adjust them at
% the end.
if mod(f*m,2) == 0
    % f*m is even; we need to cover 2*h input cells. f*m/2 should be an
    % integer, but just in case floating point roundoff causes it to
    % depart (very slight) from an integer value, snap it back with
    % round().
    h = round(f*m/2);
    c = round((a + b - 1)/2);
    first = c - (h - 1);
    last  = c + h;
else
    % f*m is odd; we need to cover 2*h + 1 input cells. As above, apply
    % round() just in case (f*m - 1)/2 takes on a slightly non-integer
    % floating point value.
    h = round((f*m - 1)/2);
    c = round((a + b)/2);
    first = c - h;
    last  = c + h;   
end

% Step 3 -- If our estimated range of input indices falls too far to the
% left, move it to the right. Likewise, if it falls too far to the
% right, move it to the left. Because we've already ensure that
% f * m <= n, only one of these conditions can hold at a time.
if first < 1
    % Shift to the right by (1 - first), which is a positive value.
    last = last + (1 - first);
    first = 1;
elseif last > n
    % Shift to the left by (last - n), which is a positive value.
    first = first - (last - n);
    last = n;
end

% Step 4 -- Re-assign the limits (in intrinsic coordinates) to match
% the extent spanned by the selected input cells.
a = first - 0.5;
b = last  + 0.5;

% Step 5 -- Now adjust first and last to correspond to the indices of the
% first and last samples to be taken from the input grid, thus matching
% their definition in the help.
if mod(f,2) == 0
    h = round(f/2);
    first = first + h - 1;
    last  = last - h;
else
    % Note: when f == 1, we end up here and h == 0.
    h = round((f - 1)/2);
    first = first + h;
    last  = last - h;
end

