function [latbin,lonbin,count] = hista(lats,lons,binarea,geoid,units)
%HISTA  Histogram for geographic points with equal-area bins
%
%  [lat,lon,ct] = HISTA(lat0,lon0) computes a spatial histogram of
%  geographic data using equal area binning.  The bin area is 100 square
%  kilometers.  The outputs are the location of bins in which the
%  data was accumulated, as well as the number of occurrences in these bins.
%
%  [lat,lon,ct] = HISTA(lat0,lon0,binarea) uses the bin size specified
%  by the input binarea, which must be in square kilometers.
%
%  [lat,lon,ct] = HISTA(lat0,lon0,binarea,geoid) assumes the data
%  is distributed on the ellipsoid defined by the input geoid.
%  The geoid vector is of the form [semimajor axes, eccentricity].
%  If omitted, the unit sphere, geoid = [1 0], is assumed.
%
%  [lat,lon,ct] = HISTA(lat0,lon0,binarea,'units') and
%  [lat,lon,ct] = HISTA(lat0,lon0,binarea,geoid,'units')
%  use the input 'units' to define the angle units of the inputs and
%  outputs.  If omitted, 'degrees' are assumed.
%
%  See also HISTR.

% Copyright 1996-2007 The MathWorks, Inc.
% $Revision: 1.13.4.8 $  $Date: 2007/11/26 20:35:15 $
% Written by:  E. Byrns, E. Brown

error(nargchk(2, 5, nargin, 'struct'))

if nargin == 2
    binarea = 100;
    geoid = [1 0];
    units = 'degrees';
elseif nargin == 3
    geoid = [1 0];
    units = 'degrees';
elseif nargin == 4
    if ischar(geoid)
        units = geoid;
        geoid = [1 0];
    else
        units = 'degrees';
    end
end

assert(isequal(size(lats),size(lons)), ...
    ['map:' mfilename ':mapError'], ...
    'Inconsistent latitude and longitude dimensions')

assert(isscalar(binarea), ...
    ['map:' mfilename ':mapError'], ...
    'Scalar bin area is required')

binarea = ignoreComplex(binarea, mfilename, 'binarea');

%  Convert to degrees and ensure column vectors
%  Ensure that the longitude data is between -180 and 180

[lats, lons] = toDegrees(units, lats(:), lons(:));
lons = npi2pi(lons,'degrees','exact');

%  Compute the mean of the input data. Center the matrix on the
%  mean to avoid problems like the north pole getting multiple 
%  bins.

datamean = meanm(lats,lons,geoid,'degrees');

%  Convert to equal area coordinates

[x,y] = grn2eqa(lats,lons,datamean,geoid,'degrees');

%  Determine the length of a side of the bin in radians

lenside = km2deg(sqrt(binarea));

%  Determine the delta in x and y direction.  Remember lenside is in radians

[xdel,ydel] = grn2eqa(lenside,lenside,[0 0 0],geoid,'degrees');

%  Determine the x and y limits of the equal area map

xlim = [min(x-xdel) max(x+xdel)];
ylim = [min(y-ydel) max(y+ydel)];

%  Construct a sparse matrix to bin the data into

[map,maplegend] = spzerom(ylim,xlim,1/xdel);

%  Bin the data into the sparse matrix
  
indx = setpostn(map,maplegend,y,x);
for i = 1:length(indx)
    map(indx(i)) = map(indx(i)) + 1;
end

%  Determine the locations of the binned data

[row,col,count] = find(map); 
[ybin,xbin]     = setltln(map,maplegend,row,col);

%  Transform the xbin and ybin back to Greenwich

[latbin,lonbin] = eqa2grn(xbin,ybin,datamean,geoid,'degrees');

%  Convert back to the proper units

[latbin, lonbin] = fromDegrees(units, latbin, lonbin);
