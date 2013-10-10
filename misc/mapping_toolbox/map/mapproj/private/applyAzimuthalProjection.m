function outputs = applyAzimuthalProjection(mproj, varargin)
%applyAzimuthalProjection  Apply an azimuthal map projection

% Copyright 2006-2010 The MathWorks, Inc.
% $Revision: 1.1.6.7 $  $Date: 2010/03/04 16:22:53 $

n = numel(varargin);
if n == 1
    outputs{1} = mproj.default(varargin{1});
elseif n == 5 || n == 6 || n == 7
    mproj.applyForward = @applyForward;
    mproj.applyInverse = @applyInverse;
    outputs = doApplyProj(mproj, varargin{:});
else
    error(['map:' mfilename ':invalidArgCount'], ...
        'Incorrect number of arguments.')
end

%--------------------------------------------------------------------------

function [x, y, savepts] = applyForward(mproj, mstruct, lat, lon, objectType)

% Convert origin to radians & auxiliary sphere
origin = convertOrigin(mproj, mstruct);

% Radius for trimming
maxRange = flatlimitToRadius(mproj, mstruct);

% Convert input coordinates to radians
[lat, lon] = toRadians(mstruct.angleunits, real(lat), real(lon));

% Convert to auxiliary latitude
lat = convertlat(mstruct.geoid, lat, ...
    'geodetic', mproj.auxiliaryLatitudeType, 'nocheck');

% Convert lat-lon to range-azimuth and trim data
newTrimming = strcmp(objectType,'notrim') || strncmp(objectType, 'geo', 3);
if newTrimming
    % Rotate before trimming
    [lat,lon] = rotatePureTriax(lat, lon, origin, 'forward');
    
    % Trim to maxRange on rotated auxiliary sphere (skip savepts)
    switch(objectType)
        
        case {'geopoint', 'geomultipoint'}
            [rng, az] = greatcircleinv00(lat, lon);
            needsTrim = (rng > maxRange);
            az( needsTrim) = [];
            rng(needsTrim) = [];
            
        case 'geoline'
            [rng, az] = greatcircleinv00(lat, lon);
            [rng, az] = trimPolylineToCircle(rng, az, maxRange);
            
        case 'geopolygon'
            inc = 2.0*pi/180;
            [lat,lon] = trimPolygonToSmallCircle(lat, lon, 0, 0, maxRange, inc);
            [rng, az] = greatcircleinv00(lat, lon);
            
        case 'geosurface'            
            % Ranges and azimuths on rotated auxiliary sphere
            [rng, az] = greatcircleinv00(lat, lon);
            
            % Trim data exceeding specified range from the origin,
            % save structure of altered points
            [rng, az] = ...
                trimdata(rng, [-inf maxRange], az,[-inf inf], 'surface');
            
        otherwise % 'notrim'
            [rng, az] = greatcircleinv00(lat, lon);
    end
    savepts.trimmed = [];
    savepts.clipped = [];
else
    % Old approach to trimming
    
    % Back off at +/- 180 degrees
    epsilon = 5*epsm('radians');
    lon = backOffAtPi(lon, epsilon);
    
    % Ranges and azimuths on rotated auxiliary sphere
    [lat,lon] = rotatemRadians(lat, lon, origin, 'forward');
    [rng, az] = greatcircleinv00(lat, lon);
    
    % Trim data exceeding specified range from the origin,
    % save structure of altered points
    [rng, az, trimmed] = ...
        trimdata(rng, [-inf maxRange], az,[-inf inf], objectType);
    savepts.trimmed = addLatLonToTrimList(trimmed, lat, lon);
    savepts.clipped = [];
end

% Project
[x, y] = mproj.forward(mstruct, rng, az);
[x, y] = applyScaleAndOriginShift(mstruct, x, y);

%--------------------------------------------------------------------------

function [lat, lon, savepts] = applyInverse( ...
    mproj, mstruct, x, y, objectType, savepts)

% Undo projection
[x, y] = undoScaleAndOriginShift(mstruct, x, y);
[rng, az] = mproj.inverse(mstruct, x, y);
[lat, lon] = reckon(0, 0, rng, az, 'radians');

% Undo trimming
[lat, lon] = undotrim(lat, lon, savepts.trimmed, objectType);

% Undo rotation and restore geodetic latitude
origin = convertOrigin(mproj, mstruct);
[lat, lon] = rotatemRadians(lat, lon, origin, 'inverse');
lat = convertlat(mstruct.geoid, lat, ...
    mproj.auxiliaryLatitudeType, 'geodetic', 'nocheck');

% Restore angle units
[lat, lon] = fromRadians(mstruct.angleunits, lat, lon);

%--------------------------------------------------------------------------

function origin = convertOrigin(mproj, mstruct)
% Convert origin to radians & auxiliary sphere

origin = toRadians(mstruct.angleunits, mstruct.origin);

origin(1) = convertlat(mstruct.geoid, origin(1), ...
    'geodetic', mproj.auxiliaryLatitudeType, 'nocheck');

%--------------------------------------------------------------------------

function radius = flatlimitToRadius(mproj, mstruct)
% Convert latitude frame limits to frame radius for azimuthal projection

% Extract the latitude frame limits and convert to radians
flatlimit = toRadians(mstruct.angleunits, mstruct.flatlimit);

% Adjust the latitudes to the auxiliary sphere
flatlimit = convertlat(mstruct.geoid, flatlimit, ...
    'geodetic', mproj.auxiliaryLatitudeType, 'nocheck');

% The maximum range (radius of circle to trim to) is encoded in flatlimit
% (usually as the second element).
radius = max(flatlimit);

%-----------------------------------------------------------------------

function [rng,az] = greatcircleinv00(phi, lambda)
% Range and azimuth from the point (0,0)

a = sin(phi/2).^2 + cos(phi) .* sin(lambda/2).^2;
rng = 2 * atan2(sqrt(a),sqrt(1 - a));
az = atan2(cos(phi) .* sin(lambda), sin(phi));
