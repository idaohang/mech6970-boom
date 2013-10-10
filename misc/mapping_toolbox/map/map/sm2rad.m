function rad=sm2rad(sm,sphere)
%SM2RAD Convert distance from statute miles to radians
%
%   RAD = SM2RAD(SM) converts distances from statute miles to radians as
%   measured along a great circle on a sphere with a radius of 6371 km
%   (3958.748 sm), the mean radius of the Earth.
%
%   RAD = SM2RAD(SM,RADIUS) converts distances from statute miles to
%   radians as measured along a great circle on a sphere having the
%   specified radius. RADIUS must be in units of statute miles.
%
%   RAD = SM2RAD(SM,SPHERE) converts distances from statute miles to
%   radians, as measured along a great circle on a sphere approximating an
%   object in the Solar System.  SPHERE may be one of the following
%   strings: 'sun', 'moon', 'mercury', 'venus', 'earth', 'mars', 'jupiter',
%   'saturn', 'uranus', 'neptune', or 'pluto', and is case-insensitive.
%
%  See also RAD2SM, SM2DEG, SM2NM, SM2KM.

% Copyright 1996-2007 The MathWorks, Inc.
% $Revision: 1.10.4.3 $  $Date: 2007/03/27 19:12:44 $

km = sm2km(sm);

if nargin == 1
    rad = km2rad(km);
else
    if ~ischar(sphere)
        sphere = sm2km(sphere);
    end
    rad = km2rad(km, sphere);
end
