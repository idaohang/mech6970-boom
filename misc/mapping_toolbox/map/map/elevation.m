function [elevationAngle, slantRange, azimuthAngle] = elevation(varargin)
%ELEVATION Local vertical elevation angle, range, and azimuth
%
%   [ELEVATIONANGLE, SLANTRANGE, AZIMUTHANGLE] = ELEVATION(LAT1, LON1, ...
%   ALT1, LAT2, LON2, ALT2) computes the elevation angle, slant range, and
%   azimuth angle of point 2 (with geodetic coordinates LAT2, LON2, and
%   ALT2) as viewed from point 1 (with geodetic coordinates LAT1, LON1, and
%   ALT1).  ALT1 and ALT2 are ellipsoidal heights.  The elevation angle is
%   the angle of the line of sight above the local horizontal at point 1.
%   The slant range is the three-dimensional Cartesian distance between
%   point 1 and point 2.  The azimuth is the angle from north to the
%   projection of the line of sight on the local horizontal. Angles are in
%   units of degrees, altitudes and distances are in meters. The figure of
%   the earth is the default ellipsoid (GRS 80) as defined by ALMANAC.
%
%   Inputs can be vectors of points, or arrays of any shape, but must match
%   in size, with the following exception:  Elevation, range, and azimuth
%   from a single point to a set of points can be computed very efficiently
%   by providing scalar coordinate inputs for point 1 and vectors or arrays
%   for point 2.
%
%   [...] = ELEVATION(LAT1,LON1, ALT1, LAT2, LON2, ALT2, ANGLEUNITS) uses
%   the string ANGLEUNITS to specify the units of the input and output
%   angles.  If omitted, 'degrees' is assumed.
%
%   [...] = ELEVATION(LAT1, LON1, ALT1, LAT2, LON2, ALT2, ANGLEUNITS,...
%   DISTANCEUNITS) uses the string DISTANCEUNITS to specify the altitude
%   and slant-range units.  If omitted, 'meters' is assumed.  Any units
%   string recognized by UNITSRATIO may be used.
% 
%   [...] = ELEVATION(LAT1, LON1, ALT1, LAT2, LON2, ALT2, ANGLEUNITS,...
%   ELLIPSOID) uses the vector ELLIPSOID, with form [semimajor axis,
%   eccentricity], to specify the ellipsoid.  If ELLIPSOID is supplied, the
%   altitudes must be in the same units as the semimajor axis and the slant
%   range will be returned in these units.  If ELLIPSOID is omitted, the
%   default earth ellipsoid defined by AZIMUTH is used and distances are in
%   meters unless otherwise specified.
%
%   Note
%   ----
%   The line-of-sight azimuth angles returned by ELEVATION will generally
%   differ slightly from the corresponding outputs of AZIMUTH and DISTANCE,
%   except for great-circle azimuths on a spherical earth. 
%
%   See also ALMANAC, AZIMUTH, DISTANCE.

% Copyright 1999-2007 The MathWorks, Inc.
% $Revision: 1.6.4.8 $  $Date: 2007/11/09 20:23:33 $

% Parse inputs, reshape input arrays into column vectors,
% convert to radians.
[phi1, lambda1, h1, phi2, lambda2, h2, angleUnits, ellipsoid, inputSize]...
    = parseInputs(varargin{:});

% Perform coordinate transformations.
[x,y,z] = geodetic2ecef(phi2,lambda2,h2,ellipsoid);
if numel(phi1) == 1
    % vectorized computation for scalar point 1
    [x,y,z] = ecef2lv(x,y,z,phi1,lambda1,h1,ellipsoid);
else
    % iteration required for non-scalar point 1
    for k = 1:numel(phi1)
       [x(k), y(k), z(k)] = ecef2lv(...
                            x(k),y(k),z(k),phi1(k),lambda1(k),h1(k),ellipsoid);
   end
end

% Convert to spherical coordinates but with azimuth defined clockwise
% from north (y) rather than counter clockwise from east (x) and in the
% interval [0 2*pi).  This the same as the computation in CART2SPH,
% except for the azimuth convention.
r = hypot(x,y);
slantRange = hypot(r,z);
elevationAngle = atan2(z,r);
azimuthAngle = mod(atan2(x,y),2*pi);

% The azimuth is undefined when point 1 is at a pole, so we choose a
% convention: zero at the north pole and pi at the south pole.
azimuthAngle(phi1 <= -pi/2) = 0;
azimuthAngle(phi1 >=  pi/2) = pi;

% Convert output angle units back from radians to match inputs.
[elevationAngle, azimuthAngle] ...
    = fromRadians(angleUnits, elevationAngle, azimuthAngle);

% Reshape outputs to match inputs.
elevationAngle = reshape(elevationAngle, inputSize);
slantRange     = reshape(slantRange,     inputSize);
azimuthAngle   = reshape(azimuthAngle,   inputSize);

%----------------------------------------------------------------------------------
function [lat1, lon1, alt1, lat2, lon2, alt2, angleUnits, ellipsoid, inputSize] ...
    = parseInputs(lat1, lon1, alt1, lat2, lon2, alt2, angleUnits, in8)

error(nargchk(6, 8, nargin, 'struct'))

% assign angle units
if nargin < 7
    angleUnits = 'degrees';
end

% assign ellipsoid/distance units
if nargin < 8
    ellipsoid = almanac('earth','geoid','m');
else
    if ischar(in8)
        % in8 is the DISTANCEUNITS string
        ellipsoid = almanac('earth','geoid',in8);
    else
        % in8 is the ELLIPSOID vector
        ellipsoid = checkellipsoid(in8,mfilename,'ELLIPSOID',8);
    end
end

% check sizes
allConsistent = isequal(...
    size(lat1),size(lon1),size(alt1),size(lat2),size(lon2),size(alt2));
scalarPoint1 = (numel(lat1) == 1 && numel(lon1) == 1 && numel(alt1) == 1);
consistentPoint2 = isequal(size(lat2), size(lon2), size(alt2));
if ~(allConsistent || (scalarPoint1 && consistentPoint2))        
    error(['map:' mfilename ':inputSizeMismatch'], ...
        'Sizes of lat, lon and alt inputs are mismatched.')
end

% convert angles to radians/coordinate arrays to column vectors
inputSize = size(lat2);
[lat1, lon1, lat2, lon2] ...
    = toRadians(angleUnits,lat1(:), lon1(:), lat2(:), lon2(:));
alt1 = alt1(:);
alt2 = alt2(:);
