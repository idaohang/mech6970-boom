function varargout = eqdcylin(varargin)
%EQDCYLIN  Equidistant Cylindrical Projection
%
%  This is a projection onto a cylinder secant at the standard parallels.
%  Distortion of both shape and area increase with distance from the
%  standard parallels.  Scale is true along all meridians (i.e. it is
%  equidistant) and the standard parallels, and is constant along any
%  parallel and along the parallel of opposite sign.
%
%  This projection was first used by Marinus of Tyre, about A.D. 100.
%  Special forms of this projection are the Plate Carree, with a standard
%  parallel at 0 deg, and the Gall Isographic, with standard parallels at
%  45 deg N and S.  Other names for this projection include
%  Equirectangular, Rectangular, Projection of Marinus, La Carte
%  Parallelogrammatique, and Die Rechteckige Plattkarte.

% Copyright 1996-2007 The MathWorks, Inc.
% $Revision: 1.9.4.4 $  $Date: 2007/03/27 19:13:49 $

mproj.default = @eqdcylinDefault;
mproj.forward = @eqdcylinFwd;
mproj.inverse = @eqdcylinInv;
mproj.auxiliaryLatitudeType = 'rectifying';
mproj.classCode = 'Cyln';

varargout = applyProjection(mproj, varargin{:});

%--------------------------------------------------------------------------

function mstruct = eqdcylinDefault(mstruct)

[mstruct.trimlat, mstruct.trimlon, mstruct.mapparallels] ...
    = fromDegrees(mstruct.angleunits, [-90 90], [-180 180], 30);
mstruct.nparallels   = 1;
mstruct.fixedorient  = [];

%--------------------------------------------------------------------------

function [x, y] = eqdcylinFwd(mstruct, lat, lon)

[radius, phi1] = deriveParameters(mstruct);

x = radius * lon * cos(phi1);
y = radius * lat;

%--------------------------------------------------------------------------

function [lat, lon] = eqdcylinInv(mstruct, x, y)

[radius, phi1] = deriveParameters(mstruct);

lat = y / radius;
lon = x / (radius*cos(phi1));

%--------------------------------------------------------------------------

function [radius, phi1] = deriveParameters(mstruct)

radius = rsphere('rectifying', mstruct.geoid(1));
phi1 = toRadians(mstruct.angleunits, mstruct.mapparallels(1));
phi1 = convertlat(...
    mstruct.geoid, phi1, 'geodetic', 'rectifying', 'nocheck');
