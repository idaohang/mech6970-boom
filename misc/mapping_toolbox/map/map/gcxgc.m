function [newlat,newlong]=gcxgc(lat1,long1,az1,lat2,long2,az2,units)
%GCXGC  Intersection points for pairs of great circles
%
%  [lat,lon] = GCXGC(lat1,long1,az1,lat2,long2,az2) finds the two
%  intersection points for every input pair of great circles.  Inputs
%  are in great circle notation, each circle defined by the latitude
%  and longitude of a point on that circle and an azimuth at that point.
%  When two circles are identical (which is not generally obvious from
%  the input notation), NaN's are returned.  Only spherical geoids are
%  supported.
%
%  [lat,lon] = GCXGC(lat1,long1,az1,lat2,long2,az2,'units') uses the
%  input string units to define the angle units for the inputs and
%  outputs.
%
%  mat = GCXGC(...) returns a single output, where mat = [lat lon].
%
%  See also SCXSC, CROSSFIX, GCXSC, RHXRH, POLYXPOLY.

% Copyright 1996-2007 The MathWorks, Inc.
% Written by:  E. Brown, E. Byrns
% $Revision: 1.11.4.4 $    $Date: 2007/11/09 20:23:59 $

error(nargchk(6, 7, nargin, 'struct'))

if nargin == 6
    units = 'degrees';
end

%  Input dimension tests

if any([length(size(lat1)) length(size(long1)) length(size(az1))  ...
        length(size(lat2)) length(size(long2)) length(size(az2))] > 2)
     error(['map:' mfilename ':mapError'], ...
         'Input matrices can not contain pages')

elseif ~isequal(size(lat1),size(long1),size(az1),...
                size(lat2),size(long2),size(az2))
	     error(['map:' mfilename ':mapError'], ...
             'Inconsistent dimensions on inputs')
end

%  Ensure real input

lat1 = real(lat1);
lat2 = real(lat2);
long1 = real(long1);
long2 = real(long2);
az1 = real(az1);
az2 = real(az2);

%  Convert input angle to radians

[lat1, lat2, long1, long2, az1, az2] ...
    = toRadians(units, lat1, lat2, long1, long2, az1, az2);

%  Convert great circle definition to a center and radius format

[lat1,long1,range1] = gc2sc(lat1,long1,az1,'radians');
[lat2,long2,range2] = gc2sc(lat2,long2,az2,'radians');

%  Compute the intersection points

[nlat,nlong] = scxsc(lat1,long1,range1,lat2,long2,range2,'radians');

%  Transform the output to the proper units

[newlat, newlong] = fromRadians(units, nlat, nlong);

%  Set the output argument if necessary

if nargout < 2;  newlat = [newlat newlong];  end
