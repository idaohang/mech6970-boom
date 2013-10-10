function km = rad2km(rad,sphere)
%RAD2KM Convert distance from radians to kilometers
%
%   KM = RAD2KM(RAD) converts distances from radians to kilometers as
%   measured along a great circle on a sphere with a radius of 6371 km, the
%   mean radius of the Earth.
%
%   KM = RAD2KM(RAD,RADIUS) converts distances from radians to kilometers
%   as measured along a great circle on a sphere having the specified
%   radius. RADIUS must be in units of kilometers.
%
%   KM = RAD2KM(RAD,SPHERE) converts distances from radians to kilometers,
%   as measured along a great circle on a sphere approximating an object in
%   the Solar System.  SPHERE may be one of the following strings: 'sun',
%   'moon', 'mercury', 'venus', 'earth', 'mars', 'jupiter', 'saturn',
%   'uranus', 'neptune', or 'pluto', and is case-insensitive.
%
%  See also KM2RAD, RADTODEG, RAD2NM, RAD2SM.

% Copyright 1996-2009 The MathWorks, Inc.
% $Revision: 1.10.4.5 $  $Date: 2009/12/28 04:17:00 $

if nargin == 1
    radius = earthRadius('km');
elseif nargin == 2 && ischar(sphere)
    radius = almanac(sphere,'radius','km');
else
    radius = sphere;
end

km = radius * rad;
