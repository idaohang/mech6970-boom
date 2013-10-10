function [latmean,lonmean] = meanm(lat,lon,geoid,units)
%MEANM  Mean location of geographic points
%
%  [latbar,lonbar] = MEANM(lat,lon) computes means for geographic
%  data.  This corresponds to the geographic centroid of a data set,
%  assuming that the data are distributed on a sphere.  In contrast,
%  MEAN assumes that the data are distributed on a cartesian plan.
%  When lat and lon are vectors, a single mean location is returned.
%  When lat and long are matrices, latmean and lonmean are row vectors
%  providing the mean locations for each column of lat and lon.
%  N-dimensional arrays are not allowed.
%
%  [latbar,lonbar] = MEANM(lat,lon,geoid) computes the geographic mean
%  on the ellipsoid defined by the input geoid. The geoid vector
%  is of the form [semimajor axes, eccentricity].  If omitted, the
%  unit sphere, geoid = [1 0], is assumed.
%
%  [latbar,lonbar] = MEANM(lat,lon,'units') use the input 'units'
%  to define the angle units of the inputs and outputs.  If
%  omitted, 'degrees' are assumed.
%
%  [latbar,lonbar] = MEANM(lat,lon,geoid,'units') is a valid form.
%
%  mat = MEANM(...) returns a single output, where mat = [latbar,lonbar].
%  This is particularly useful if the lat and lon inputs are vectors.
%
%  See also MEAN, STDM, STDIST.

% Copyright 1996-2007 The MathWorks, Inc.
% $Revision: 1.9.4.5 $  $Date: 2007/11/09 20:24:41 $
% Written by:  E. Byrns, E. Brown, W. Stumpf

error(nargchk(2, 4, nargin, 'struct'))

if nargin == 2
    geoid = [];       units = [];
elseif nargin == 3
    if ischar(geoid)
        units = geoid;
        geoid = [];
    else
        units = [];
    end
end

%  Empty argument tests

if isempty(units);   units = 'degrees';   end
if isempty(geoid);   geoid = [1 0];       end

%  Input dimension tests

if ndims(lat)>2
	error(['map:' mfilename ':mapError'], ...
        'Latitude and longitude inputs limited to two dimensions.')
end

if ~isequal(size(lat),size(lon))
    error(['map:' mfilename ':mapError'], ...
        'Inconsistent dimensions on latitude and longitude inputs')
end

%  Test the geoid parameter

geoid = geoidtst(geoid);

%  Convert inputs to radians.  Use an authalic sphere for
%  calculations.  Thus, the mean calculation is area-based,
%  since the authalic sphere has the same surface area as
%  the ellipsoid.

[lat, lon] = toRadians(units, lat, lon);
lat = convertlat(geoid, lat, 'geodetic', 'authalic', 'nocheck');

%  Convert the input data to cartesian coordinates.
%  Compute the centroid point by vector summing all cartesian data.

[x,y,z]=sph2cart(lon,lat,ones(size(lat)));
[lonbar,latbar,radius]=cart2sph(sum(x),sum(y),sum(z));

%  Transform outputs to proper units.  Set longitude in -pi to pi range

latbar = convertlat(geoid, latbar, 'authalic', 'geodetic', 'nocheck');
lonbar = npi2pi(lonbar,'radians','exact');

[latbar, lonbar] = fromRadians(units, latbar, lonbar);

%  Eliminate any points whose vector sum produces a point
%  which is near the center of the sphere or ellipsoid.  This
%  occurs when the data consists of only points and their
%  antipodes, and in this case, the centroid does not have any meaning.

indx = find(radius <= epsm('radians'));
if ~isempty(indx);
    warning('map:meanm:allPointsCancel', ...
'Data in at least one column consists of only points and their antipodes.')
	latbar(indx) = NaN;   lonbar(indx) = NaN;
end

%  Set the output arguments

if nargout < 2
    latmean = [latbar lonbar];
elseif nargout == 2
    latmean = latbar;   lonmean = lonbar;
end
