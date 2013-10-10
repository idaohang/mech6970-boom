function dist=departure(long1,long2,lat,geoid,units)
%DEPARTURE  Departure of longitudes at specific latitudes
%
%  d = DEPARTURE(long1,long2,lat) computes the departure distance from
%  long1 to long2 at the input latitude lat.  Departure is the
%  distance along a specific parallel between two meridians.  The
%  output d is returned in degrees of arc length on a sphere.
%
%  d = DEPARTURE(long1,long2,lat,geoid) computes the departure
%  assuming that the input points lie on the ellipsoid defined by
%  the input geoid.  The geoid vector is of the form
%  [semimajor axes, eccentricity].
%
%  d = DEPARTURE(long1,long2,lat,'units') uses the input string 'units'
%  to define the angle units of the input and output data.  In this
%  form, the departure is returned as an arc length in the units
%  specified by 'units'.  If 'units' is omitted, 'degrees' are assumed.
%
%  d = DEPARTURE(long1,long2,lat,geoid,'units') is a valid calling
%  form.  In this case, the departure is computed in the same units as
%  the semimajor axes of the geoid vector.
%
%  See also DISTANCE.

% Copyright 1996-2007 The MathWorks, Inc.
% $Revision: 1.9.4.5 $  $Date: 2007/11/26 20:35:11 $
% Written by:  E. Byrns, E. Brown

error(nargchk(3, 5, nargin, 'struct'))

if nargin == 3
    geoid = [1 0];
    nogeoid = 1;
    units = 'degrees';
elseif nargin == 4
    if ischar(geoid)
        units = geoid;
        geoid = [1 0];
        nogeoid = 1;
    else
        units = 'degrees';
        nogeoid = 0;
    end
else
    nogeoid = 0;
end

%  Input dimension tests

if ~isequal(size(long1),size(long2))
    error(['map:' mfilename ':mapError'], ...
        'Inconsistent dimensions on longitude inputs')
end

if ~isequal(size(long1),size(lat))
    error(['map:' mfilename ':mapError'], ...
        'Inconsistent dimensions on latitude and longitude inputs')
end

%  Test the geoid parameter
geoid = geoidtst(geoid);

%  Convert inputs to radians.  Ensure that longitudes are
%  in the range 0 to 2pi since they will be treated as distances.

[lat, long1, long2] = toRadians(units, lat, long1, long2);

long1 = zero22pi(long1,'radians','exact');
long2 = zero22pi(long2,'radians','exact');

r=rcurve('parallel',geoid,lat,'radians');

deltalong=abs(long1-long2);

dist=r.*deltalong;


if nogeoid
	dist = fromRadians(units, dist);
end
