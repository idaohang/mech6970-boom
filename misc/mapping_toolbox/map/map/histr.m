function [latbin,lonbin,count,wcount] = histr(lats,lons,bindensty,units)
%HISTR  Histogram for geographic points with equirectangular bins
%
%  [lat,lon,ct] = HISTR(lat0,lon0) computes a spatial histogram of
%  geographic data using equirectangular binning of one degree.  In
%  other words, one degree increments of latitude and longitude to
%  define the bins throughout the globe.  As a result, these bins are
%  not equal area.  The outputs are the location of bins in which the
%  data was accumulated, as well as the number of occurrences in these bins.
%
%  [lat,lon,ct] = HISTR(lat0,lon0,binsize) uses the bin size specified
%  by the input.  This input must be in the same units as the lat and
%  lon input.
%
%  [lat,lon,ct] = HISTR(lat0,lon0,'units') and
%  [lat,lon,ct] = HISTR(lat0,lon0,binsize,'units') use the input
%  'units' to define the angle units of the inputs and outputs.
%  If omitted, 'degrees' are assumed.
%
%  [lat,lon,ct,wt] = HISTR(...) returns the number of occurrences,
%  weighted by the area of each bin.  The weighting factors assume that
%  bins along the equator are given an area of 1.0.
%
%  See also HISTA.

% Copyright 1996-2007 The MathWorks, Inc.
% $Revision: 1.10.4.7 $  $Date: 2007/11/26 20:35:18 $
% Written by:  E. Byrns, E. Brown

error(nargchk(2, 4, nargin, 'struct'))

if nargin == 2
    units = 'degrees';
    bindensty = 1;
elseif nargin == 3
    if ischar(bindensty)
        units = bindensty;
        bindensty = 1;
    else
        units = 'degrees';
    end
end

assert(isequal(size(lats),size(lons)), ...
    ['map:' mfilename ':mapError'], ...
    'Inconsistent latitude and longitude dimensions')

assert(isscalar(bindensty), ...
    ['map:' mfilename ':mapError'], ...
    'Scalar bin density is required')

bindensty = ignoreComplex(bindensty, mfilename, 'bindensty');

%  Convert to degrees and ensure column vectors
%  Ensure that the longitude data is between -180 and 180

[lats, lons] = toDegrees(units, lats(:), lons(:));
lons = wrapTo180(lons);

%  Construct a sparse matrix to bin the data into

latlim = [floor(min(lats))   ceil(max(lats))];
lonlim = [floor(min(lons))   ceil(max(lons))];

[map,maplegend] = spzerom(latlim,lonlim,bindensty);

%  Bin the data into the sparse matrix

indx = setpostn(map,maplegend,lats,lons);
for i = 1:length(indx)
    map(indx(i)) = map(indx(i)) + 1;
end

%  Determine the locations of the binned data

[row,col,count] = find(map);
[latbin,lonbin] = setltln(map,maplegend,row,col);

%  Convert back to the proper units

[latbin, lonbin] = fromDegrees(units, latbin, lonbin);

%  Determine the data occurrences weighted by the bin area
%  If this output is not requested, don't waste time calculating it.

if nargout == 4
    [a,areavec]=areamat(map>0,maplegend);
    wcount = full(max(areavec(row)) * count ./ areavec(row) );
end
