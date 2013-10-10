function a=areaquad(lat1,lon1,lat2,lon2,in5,in6)
%AREAQUAD Surface area of latitude-longitude quadrangle
%
%   a = AREAQUAD(lat1,lon1,lat2,lon2) returns the surface area of the
%   geographic quadrangle bounded by the parallels lat1 and lat2 and the
%   meridians lon1 and lon2. The output area is a fraction of the unit
%   sphere's area of 4*pi, so the result ranges from 0 to 1.
%
%   a = AREAQUAD(lat1,lon1,lat2,lon2,ELLIPSOID) allows the specification of
%   the ellipsoid model with the two-element ellipsoid vector ELLIPSOID.
%   When ELLIPSOID is input, the resulting area is given in terms of the
%   square of the unit of length used to define the ellipsoid. For example,
%   if the ellipsoid ALMANAC('earth','ellipsoid','kilometers') is used, the
%   resulting area is in square kilometers. The default ellipsoid is the
%   unit sphere.
%
%   a = AREAQUAD(lat1,lon1,lat2,lon2,ELLIPSOID,'angleunits') specifies the
%   angle units of the inputs. The default is 'degrees'.
%
%  See also AREAINT, AREAMAT.

% Copyright 1996-2010 The MathWorks, Inc.
% Written by:  E. Brown, E. Byrns
% $Revision: 1.9.4.6 $  $Date: 2010/09/24 14:33:07 $

error(nargchk(4,6,nargin,'struct'))

if nargin==4
	units = [];
    ellipsoid = [];
elseif nargin==5
	if ischar(in5)
		units = in5;
        ellipsoid = [];
	else
		units = [];
        ellipsoid = in5;
	end
elseif nargin==6
	ellipsoid=in5;
    units=in6;
end

%  Empty argument tests

if isempty(units)
    units  = 'degrees';
end

absolute_units = 1;          %  Report answer in absolute units assuming
if isempty(ellipsoid)        %  a radius input has been supplied.  Otherwise,
     ellipsoid = [1 0];      %  report surface area answer as a fraction
	 absolute_units = 0;     %  of a sphere
end

%  Test the ellipsoid parameter
ellipsoid = checkellipsoid(ellipsoid,mfilename,'ELLIPSOID',5);

%  Input dimension tests
validateattributes(lat1,{'double'},{'real'},mfilename,'LAT1',1);
validateattributes(lon1,{'double'},{'real'},mfilename,'LON1',2);
validateattributes(lat2,{'double'},{'real'},mfilename,'LAT2',3);
validateattributes(lon2,{'double'},{'real'},mfilename,'LON2',4);

if ~isequal(size(lat1),size(lon1),size(lat2),size(lon2))
    error('map:areaquad:latlonSizeMismatch', ...
        'Latitude and longitude inputs must all match in size.');
end

%  Convert angles to radians and transform to the authalic sphere

[lat1,lon1,lat2,lon2] = toRadians(units,lat1,lon1,lat2,lon2);
lat1 = convertlat(ellipsoid, lat1, 'geodetic', 'authalic', 'nocheck');
lat2 = convertlat(ellipsoid, lat2, 'geodetic', 'authalic', 'nocheck');
radius = rsphere('authalic',ellipsoid);

%  Compute the surface area as a fraction of a sphere

a = abs(lon1-lon2) .* abs(sin(lat1)-sin(lat2)) / (4*pi);

%  Convert to absolute terms if the default radius was not used

if absolute_units;
    a = a * 4*pi*radius^2;
end
