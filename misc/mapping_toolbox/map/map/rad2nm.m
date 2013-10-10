function nm = rad2nm(rad,sphere)
%RAD2NM Convert distance from radians to nautical miles
%
%   NM = RAD2NM(RAD) converts distances from radians to nautical miles as
%   measured along a great circle on a sphere with a radius of 6371 km 
%   (3440.065 nm), the mean radius of the Earth.
%
%   NM = RAD2NM(RAD,RADIUS) converts distances from radians to nautical
%   miles as measured along a great circle on a sphere having the specified
%   radius. RADIUS must be in units of nautical miles.
%
%   NM = RAD2NM(RAD,SPHERE) converts distances from radians to nautical
%   miles, as measured along a great circle on a sphere approximating an
%   object in the Solar System.  SPHERE may be one of the following
%   strings: 'sun', 'moon', 'mercury', 'venus', 'earth', 'mars', 'jupiter',
%   'saturn', 'uranus', 'neptune', or 'pluto', and is case-insensitive.
%
%  See also NM2RAD, RADTODEG, RAD2KM, RAD2SM.

% Copyright 1996-2009 The MathWorks, Inc.
% $Revision: 1.10.4.4 $  $Date: 2009/03/30 23:38:43 $

if nargin == 1
    km = rad2km(rad);
else
    if ~ischar(sphere)
        sphere = nm2km(sphere);
    end
    km = rad2km(rad,sphere);
end

nm = km2nm(km);
