function [az,len] = vinvtran(varargin)
%VINVTRAN  Azimuth on ellipsoid from direction angle in map plane
%
%  az = VINVTRAN(x,y,th) transforms an angle in the projection space at
%  the point specified by x and y  into an azimuth angle in Greenwich
%  coordinates. The map projection currently displayed is used to 
%  define the projection space.  The input angles must be in the same 
%  units as specified by the current map projection.  The inputs can be 
%  scalars or matrices of the equal size. The angle in the projection space
%  angle th is defined positive counter-clockwise from the x axis.
%
%  az = VINVTRAN(struct,x,y,th) uses the map projection defined by the
%  input struct to compute the map projection.
% 
%  [az,len] = VINVTRAN(...) also returns the vector length in the 
%  Greenwich coordinate system. A value of 1 indicates no scale distortion
%  for that angle.
%
%  This transformation is limited to the region specified by
%  the frame limits in the current map definition.
%
%  See also VFWDTRAN, MFWDTRAN, MINVTRAN, DEFAULTM.

% Copyright 1996-2008 The MathWorks, Inc.
% $Revision: 1.1.6.5 $  $Date: 2008/02/07 16:31:03 $
% Written by:  E. Byrns, W. Stumpf

error(nargchk(3, 4, nargin, 'struct'))

%  Parse inputs

if nargin == 3
    mstruct = gcm;
else
    % nargin == 4
    if isstruct(varargin{1})
        mstruct = varargin{1};
        varargin(1) = [];
    else
        eid = sprintf('%s:%s:missingstruct', getcomp, mfilename);
        error(eid, ...
            'When four arguments are supplied, the first must be\na map projection structure.)')
    end
end

x   = real(varargin{1});
y   = real(varargin{2});
th  = real(varargin{3});

%  Check inputs

if ~isequal(size(x),size(y),size(th))
    eid = sprintf('%s:%s:inconsistentInputSizes', getcomp, mfilename);
    error(eid, 'x, y and theta inputs must be the same size.')
end

if strcmp(mstruct.mapprojection,'globe')
    eid = sprintf('%s:%s:globe', getcomp, mfilename);
    error(eid, 'VINVTRAN does not work on globe projections.')
end

epsilon = 100*epsm('degrees');
units = mstruct.angleunits;

%  Transform x and y to lat and long

[lat,lon] = minvtran(mstruct,x,y);

%  Transform data to degrees

[lat, lon, origin, frmlon, frmlat] = toDegrees(units, ...
    lat, lon, mstruct.origin, mstruct.flonlimit, mstruct.flatlimit);

%  Rotate the input data to the base coordinate system.
%  This is the same coordinate system as the map frame.

[LatRot,LonRot] = rotatem(lat,lon,origin,'forward','degrees');

%  Check for points outside the map frame

indx = (LonRot < min(frmlon) | LonRot > max(frmlon) | ...
        LatRot < min(frmlat) | LatRot > max(frmlat) );
if any(indx(:))
    wid = sprintf('%s:%s:pointOutsideFrame', getcomp, mfilename);
    warning(wid,'Point outside of valid projection region.')
    LatRot(indx) = NaN;
    LonRot(indx) = NaN;
end

%  Check for points near the edge of the map

% Back away from the poles to avoid problems reckoning. Convergence of
% the meridians makes east-west movements cross the dateline in longitude, 
% even if we back away by a couple of meters.

latlim = 89.9;
dlat = 90 - latlim;

LatRot(LatRot <= -latlim) = LatRot(LatRot <= -latlim) + dlat;
LatRot(LatRot >=  latlim) = LatRot(LatRot >=  latlim) - dlat; 

% Back away from the edges

LonRot(LonRot <= min(frmlon) + epsilon) = min(frmlon) + epsilon;
LonRot(LonRot >= max(frmlon) - epsilon) = max(frmlon) - epsilon;
LatRot(LatRot <= min(frmlat) + epsilon) = min(frmlat) + epsilon;
LatRot(LatRot >= max(frmlat) - epsilon) = max(frmlat) - epsilon;

%  Return processed data back to the original space

[LatNew,LonNew] = rotatem(LatRot,LonRot,origin,'inverse','degrees');

%  Transform input data to the projection space

[x,y] = mfwdtran(mstruct,LatNew,LonNew);

%  Compute a new point off the starting point

%  We need to define the size of a small step, RNG, that we'll take in
%  the direction of angle TH.  On the Earth, a step size of about 10
%  meters seems reasonable.  Because it's about 10,000,000 meters
%  from Equator to either pole, we can generalize our step size to be
%  1/1,000,000 of the Equator-to-pole distance.  And given a
%  slightly-flattened ellipsoid with semi-major axis length a, the
%  Equator-to-pole distance is roughly a * pi/2.  Therefore we'll use
%  RNG = 1e-6 * a * pi/2:
rng = 1e-6 * mstruct.geoid(1) * pi/2;
th = toRadians(units, th);

x2 = x + rng .* cos(th);
y2 = y + rng .* sin(th);

%  Compute the second point in Greenwich coordinates

[Lat2,Lon2] = minvtran(mstruct,x2,y2);

%  Compute the azimuth and length of vector

[dist, az] = distance('rh',LatNew,LonNew,Lat2,Lon2,mstruct.geoid,'degrees');

%  Transform back to output units

az = fromDegrees(units, az);
len = dist ./ rng;
