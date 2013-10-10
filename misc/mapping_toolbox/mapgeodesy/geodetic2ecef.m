function [x, y, z] = geodetic2ecef(phi, lambda, h, ellipsoid)
%GEODETIC2ECEF Convert geodetic to geocentric (ECEF) coordinates
%
%   [X, Y, Z] = GEODETIC2ECEF(PHI, LAMBDA, H, ELLIPSOID) converts geodetic
%   point locations specified by the coordinate arrays PHI (geodetic
%   latitude in radians), LAMBDA (longitude in radians), and H (ellipsoidal
%   height) to geocentric Cartesian coordinates X, Y, and Z.  The geodetic
%   coordinates refer to the reference ellipsoid specified by ELLIPSOID (a
%   row vector with the form [semimajor axis, eccentricity]).  H must use
%   the same units as the semimajor axis;  X, Y, and Z will be expressed in
%   these units also.
%
%   The geocentric Cartesian coordinate system is fixed with respect to the
%   Earth, with its origin at the center of the ellipsoid and its X-, Y-,
%   and Z-axes intersecting the surface at the following points:
%
%                PHI  LAMBDA
%      X-axis:    0     0      (Equator at the Prime Meridian)
%      Y-axis:    0   pi/2     (Equator at 90-degrees East)
%      Z-axis:  pi/2    0      (North Pole)
%
%   A common synonym is Earth-Centered, Earth-Fixed coordinates, or ECEF.
%
%   See also ECEF2GEODETIC, ECEF2LV, GEODETIC2GEOCENTRICLAT, LV2ECEF.

% Copyright 2004-2009 The MathWorks, Inc.
% $Revision: 1.1.6.4 $  $Date: 2009/04/15 23:34:46 $

% Reference
% ---------
% Paul R. Wolf and Bon A. Dewitt, "Elements of Photogrammetry with
% Applications in GIS," 3rd Ed., McGraw-Hill, 2000 (Appendix F-3).

a  = ellipsoid(1);
e2 = ellipsoid(2) ^ 2;
sinphi = sin(phi);
cosphi = cos(phi);
N  = a ./ sqrt(1 - e2 * sinphi.^2);
x = (N + h) .* cosphi .* cos(lambda);
y = (N + h) .* cosphi .* sin(lambda);
z = (N*(1 - e2) + h) .* sinphi;
