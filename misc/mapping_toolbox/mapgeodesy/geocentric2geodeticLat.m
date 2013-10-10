function phi = geocentric2geodeticLat(ecc, phi_g)
%GEOCENTRIC2GEODETICLAT  Convert geocentric to geodetic latitude
%
%   PHI = GEOCENTRIC2GEODETICLAT(ECC, PHI_G) converts an array of
%   geocentric latitude in radians, PHI_G, to geodetic latitude in
%   radians, PHI, on a reference ellipsoid with first eccentricity ECC.
%   (ECC is the second element of the reference ellipsoid vector.)
%
%   For conversion to/from other types of auxiliary latitude and,
%   optionally, to work in degrees, use Mapping Toolbox function
%   CONVERTLAT.  For conversion from 3-D geocentric coordinates, see
%   ECEF2GEODETIC.
%
%   See also CONVERTLAT, ECEF2GEODETIC, GEODETIC2GEOCENTRICLAT.

% Copyright 2006 The MathWorks, Inc.
% $Revision: 1.1.8.1 $  $Date: 2006/05/17 21:04:33 $

% Technical reference:
%    J. P. Snyder, "Map Projections - A Working Manual,"  US Geological
%    Survey Professional Paper 1395, US Government Printing Office,
%    Washington, DC, 1987, pp. 13-18.

phi = atan2( sin(phi_g), (1-ecc^2)*cos(phi_g) );
