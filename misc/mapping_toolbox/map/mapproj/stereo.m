function varargout = stereo(varargin)
%STEREO  Stereographic Azimuthal Projection
%
%  This is a perspective projection on a plane tangent at the center
%  point from the point antipodal to the center point.  The center
%  point is a pole in the common polar aspect, but it can be any
%  point.  This projection has two significant properties.  It is
%  conformal, being free from angular distortion.  Additionally, all
%  great and small circles are either straight lines or circular
%  arcs on this projection.  Scale is true only at the center point,
%  and is constant along any circle having the center point as its
%  center.  This projection is not equal area.
%
%  The polar aspect of this projection appears to have been developed
%  by the Egyptians and Greeks by the second century B.C.

% Copyright 1996-2007 The MathWorks, Inc.
% $Revision: 1.10.4.4 $  $Date: 2007/03/27 19:14:18 $

mproj.default = @stereoDefault;
mproj.forward = @stereoFwd;
mproj.inverse = @stereoInv;
mproj.auxiliaryLatitudeType = 'conformal';
mproj.classCode = 'Azim';

varargout = applyAzimuthalProjection(mproj, varargin{:});

%--------------------------------------------------------------------------

function mstruct = stereoDefault(mstruct)

[mstruct.trimlat, mstruct.trimlon] ...
    = fromDegrees(mstruct.angleunits, [-Inf 90], [-180 180]);
mstruct.mapparallels = [];
mstruct.nparallels   = 0;
mstruct.fixedorient  = [];

%--------------------------------------------------------------------------

function [x, y] = stereoFwd(mstruct, rng, az)

fact1 = deriveParameters(mstruct);
r = (fact1 * sin(rng)) ./ ( 1 + cos(rng) );

x = r .* sin(az);
y = r .* cos(az);

%--------------------------------------------------------------------------

function [rng, az] = stereoInv(mstruct, x, y)

fact1 = deriveParameters(mstruct);
rho = (x.^2 + y.^2) / fact1^2;

az = atan2(x, y);
rng = acos((1 - rho) ./ (1 + rho));

%--------------------------------------------------------------------------

function fact1 = deriveParameters(mstruct)

% Derive projection parameters.

a = mstruct.geoid(1);
e = mstruct.geoid(2);

phi0 = toRadians(mstruct.angleunits, mstruct.origin(1));
chi0 = convertlat(mstruct.geoid, phi0, 'geodetic', 'conformal', 'nocheck');

den1 = (1 - (e*sin(phi0))^2);
m1 = cos(phi0) / sqrt(den1);
epsilon = epsm('radians');
if abs(pi/2 - abs(phi0)) <= epsilon
    fact1 = 2*a / sqrt(den1);
else
    fact1 = 2*a*m1 / cos(chi0);
end
