function phi_g = geodetic2geocentricLat(ecc, phi)
%GEODETIC2GEOCENTRICLAT  Convert geodetic to geocentric latitude
%
%   PHI_G = GEODETIC2GEOCENTRICLAT(ECC, PHI) converts an array of
%   geodetic latitude in radians, PHI, to geocentric latitude in
%   radians, PHI_G, on a reference ellipsoid with first eccentricity ECC
%   (ECC is the second element of the reference ellipsoid vector.)
%
%   For conversion to/from other types of auxiliary latitude and,
%   optionally, to work in degrees, use Mapping Toolbox function
%   CONVERTLAT.  For conversion to 3-D geocentric coordinates, see
%   GEODETIC2ECEF.
%
%   See also CONVERTLAT, GEOCENTRIC2GEODETICLAT, GEODETIC2ECEF.

% Copyright 2006 The MathWorks, Inc.
% $Revision: 1.1.8.1 $  $Date: 2006/05/17 21:04:34 $

% Technical reference:
%    J. P. Snyder, "Map Projections - A Working Manual,"  US Geological
%    Survey Professional Paper 1395, US Government Printing Office,
%    Washington, DC, 1987, pp. 13-18.

phi_g = atan2( (1-ecc^2)*sin(phi), cos(phi) );
