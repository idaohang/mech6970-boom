function [newlat,newlong]=gcxsc(gclat,gclong,gcaz,sclat,sclong,scrange,units)
%GCXSC  Intersection points for great and small circle pairs
%
%  [lat,lon] = GCXSC(gclat,gclong,gcaz,sclat,sclong,scrange) finds the
%  intersection points, if any, between a great circle given in great
%  circle notation (lat,long,azimuth, lat/long on the circle) and a
%  small circle given in small circle notation (lat,long,range,
%  lat/long the center of the circle). GCXSC returns NaNs if the
%  circles do not intersect or are identical.  The great circle(s)
%  must be input first.  Only spherical geoids are supported.
%
%  [lat,lon] = GCXSC(gclat,gclong,gcaz,sclat,sclong,scrange,'units')
%  uses the input string units to define the angle units for the
%  inputs and outputs.
%
%  mat = GCXSC(...) returns a single output, where mat = [lat lon].
%
%  See also SCXSC, CROSSFIX, GCXGC, RHXRH.

% Copyright 1996-2007 The MathWorks, Inc.
% $Revision: 1.10.4.4 $  $Date: 2007/11/09 20:24:00 $

error(nargchk(6, 7, nargin, 'struct'))

if nargin == 6
    units = 'degrees';
end

%  Input dimension tests

if any([ndims(gclat) ndims(gclong) ndims(gcaz)  ...
        ndims(sclat) ndims(sclong) ndims(scrange)] > 2)
     error(['map:' mfilename ':mapError'], ...
         'Input matrices can not contain pages')

elseif ~isequal(size(gclat),size(gclong),size(gcaz),...
                size(sclat),size(sclong),size(scrange))
	     error(['map:' mfilename ':mapError'], ...
             'Inconsistent dimensions on inputs')
end

%  Ensure real input

gclat = real(gclat);
sclat = real(sclat);
gclong = real(gclong);
sclong = real(sclong);
gcaz = real(gcaz);
scrange = real(scrange);

%  Convert input angle to radians

[gclat, sclat, gclong, sclong, gcaz, scrange] ...
    = toRadians(units, gclat, sclat, gclong, sclong, gcaz, scrange);

%  Convert great circle definition to a center and radius format

[gclat,gclong,gcrange]=gc2sc(gclat,gclong,gcaz,'radians');

%  Compute the intersection points

[newlat,newlong]=scxsc(gclat,gclong,gcrange,sclat,sclong,scrange,'radians');

%  Transform the output to the proper units

[newlat, newlong] = fromRadians(units, newlat, newlong);

%  Set the output argument if necessary

if nargout < 2;  newlat = [newlat newlong];  end
