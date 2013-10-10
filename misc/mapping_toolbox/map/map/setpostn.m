function varargout = setpostn(Z, R, lat, lon)
%SETPOSTN  Convert latitude-longitude to data grid rows and columns
%
%   [ROW, COL] = SETPOSTN(Z, R, LAT, LON) returns the row and column
%   indices of the regular data grid Z for the points specified by the
%   vectors LAT and LON. R can be a spatialref.GeoRasterReference
%   object, a referencing vector, or a referencing matrix.
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
%   along a meridian and each row falls along a parallel. Points falling
%   outside the grid are ignored in ROW and COL.  All input angles are
%   in degrees.
%
%   INDX = SETPOSTN(...) returns the indices of Z corresponding to
%   the points in LAT and LON.  Points falling outside the grid are
%   ignored in INDX.
%
%   [ROW, COL, indxPointOutsideGrid] = SETPOSTN(...) returns the indices
%   of LAT and LON corresponding to points outside the grid.  These
%   points are ignored in ROW and COL.
%
%   See also LATLON2PIX, LTLN2VAL, SETLTLN.

% Copyright 1996-2010 The MathWorks, Inc.
% $Revision: 1.12.4.11 $  $Date: 2010/11/17 11:24:37 $

lat = ignoreComplex(lat, mfilename, 'LAT');
lon = ignoreComplex(lon, mfilename, 'LON');
checklatlon(lat, lon, mfilename, 'LAT', 'LON', 3, 4)

%  If R is already spatial referencing object, validate it. Otherwise
%  convert the input referencing vector or matrix.
R = internal.map.convertToGeoRasterRef( ...
    R, size(Z), 'degrees', 'SETPOSTN', 'R', 2);

if strcmp(R.RasterInterpretation, 'postings')
    [row, col] = R.geographicToSub(lat, lon);
    indxPointOutsideGrid = find(isnan(row));
    row(indxPointOutsideGrid) = [];
    col(indxPointOutsideGrid) = [];
else
    %  Identify points located outside the grid limits
    indxPointOutsideGrid = find(~R.contains(lat, lon));
    
    %  Remove such points from further processing and warn
    if ~isempty(indxPointOutsideGrid)
        lat(indxPointOutsideGrid) = [];
        lon(indxPointOutsideGrid) = [];
        warning('map:setpostn:pointOutsideLimits', ...
            'At least one point falls outside of the limits of the data grid.')
    end
    
    %  Compute scaled offsets relative to minimum grid limits
    row = R.latitudeToIntrinsicY(lat) - R.latitudeToIntrinsicY(R.FirstCornerLat);
    col = R.longitudeToIntrinsicX(lon) - R.longitudeToIntrinsicX(R.FirstCornerLon);
    
    %  Convert scaled subscripts into integer row and column subscripts
    row = ceil(row);
    col = ceil(col);
    
    % Even though we've already used INGEOQUAD to eliminate all inputs that
    % fall outside the grid, it's still posssible for roundoff effects in
    % the preceding four lines to cause out-of-range results, so we clamp
    % row to [1 ... size(Z,1)] and col to [1 ... size(Z,2)].
    row = max(row,1);
    col = max(col,1);
    row = min(row,size(Z,1));
    col = min(col,size(Z,2));
end

%  Assign outputs
switch(nargout)
    case {0,1}
        varargout = {(col-1)*size(Z,1) + row};
    case 2
        varargout = {row, col};
    case 3
        varargout = {row, col, indxPointOutsideGrid};
end
