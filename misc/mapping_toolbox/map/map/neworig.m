function [mapout,lat,lon] = neworig(map,R,origin,direction,units)
%NEWORIG Orient regular data grid to oblique aspect
%
%   [Z,lat,lon] = NEWORIG(Z0,R,origin) and
%   [Z,lat,lon] = NEWORIG(Z0,R,origin,'forward') will transform regular
%   data grid Z0 into an oblique aspect, while preserving the matrix
%   storage format.  In other words, the oblique map origin is not
%   necessarily at (0,0) in the Greenwich coordinate frame. This allows
%   operations to be performed on the matrix representing the oblique map.
%   For example, azimuthal calculations for a point in a data grid become
%   row and column operations if the data grid is transformed so that the
%   north pole of the oblique map represents the desired point on the
%   globe.  R can be a spatialref.GeoRasterReference object, a
%   referencing vector, or a referencing matrix.
%
%   If R is a spatialref.GeoRasterReference object, its RasterSize
%   property must be consistent with size(Z).
%
%   If R is a referencing vector, it must be a 1-by-3 with elements:
%
%     [cells/degree northern_latitude_limit western_longitude_limit]
%
%   If R is a referencing matrix, it must be 3-by-2 and transform raster
%   row and column indices to/from geographic coordinates according to:
% 
%                     [lon lat] = [row col 1] * R.
%
%   If R is a referencing matrix, it must define a (non-rotational,
%   non-skewed) relationship in which each column of the data grid falls
%   along a meridian and each row falls along a parallel. 
%
%   [Z,lat,lon] = NEWORIG(Z0,R,origin,'inverse') transforms
%   the regular data grid from the oblique frame to the Greenwich
%   coordinate frame.
%
%   See also ROTATEM, SETPOSTN, ORG2POL.

% Copyright 1996-2010 The MathWorks, Inc.
% $Revision: 1.11.4.8 $  $Date: 2010/11/17 11:24:31 $
% Written by:  E. Brown, E. Byrns

error(nargchk(3, 5, nargin, 'struct'))

if nargin == 3
    direction = [];
    units = [];
elseif nargin == 4
    units = [];
end

%  Empty argument tests

if isempty(direction)
    direction = 'forward';
else
    direction = validatestring(direction, {'forward','inverse'}, ...
    'neworig', 'DIRECTION', 4);
end

if isempty(units) 
    units = 'degrees'; 
end

%   Compute the starting grid locations

[lat,lon] = meshgrat(map,R,size(map));

%  Convert units

[lat, lon, origin] = toRadians(units, lat, lon, origin);

%  Set the proper direction for rotatem.  If the user has entered
%  forward, then this is actually an inverse using rotatem.  We
%  must find out what Greenwich coordinates will produce the [lat lon]
%  grid.  It is the codes of these Greenwich coordinates that we
%  want to move into the positions corresponding to [lat lon].  This
%  process works in reverse when the user is going in the inverse direction.

if strcmp(direction,'forward')
     direction = 'inverse';
else
     direction = 'forward';
end

%  Rotate the grid to the corresponding starting locations.

[lat,lon] = rotatem(lat,lon,origin,direction);

%  Convert the grid to the units in which the map is stored

[lat, lon] = fromRadians(units, lat, lon);

%  Compute the starting positions of the map coordinates

indx = setpostn(map,R,lat(:),lon(:));

%  Set these indices in the new map.

[r,c] = size(map);
mapout = map(indx);
mapout = reshape(mapout,r,c);
