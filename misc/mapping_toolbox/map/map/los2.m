function [vis,visprofile,dist,z,lattrk,lontrk] ...
  = los2(map,R,lat1,lon1,lat2,lon2,oalt,talt,...
         oaltopt,taltopt,actualradius,apparentradius)
%LOS2 Line of sight visibility between two points in terrain
%
%   LOS2 computes the mutual visibility between two points on a displayed
%   digital elevation map.  LOS2 uses the current object if it is a regular
%   data grid, or the first regular data grid found on the current axes.
%   The grid's zdata is used for the profile.  The color data is used in
%   the absence of data in z.  The two points are selected by clicking on
%   the map. The result is displayed in a new figure.  Markers indicate
%   visible and obscured points along the profile.  The profile is shown in
%   a Cartesian coordinate system with the origin at the observer's
%   location.  The displayed z coordinate accounts for the elevation of the
%   terrain and the curvature of the body.
%
%   VIS = LOS2(Z, R, LAT1, LON1, LAT2, LON2) computes the mutual visibility
%   between pairs of points on a digital elevation map.  The elevations are
%   provided as a regular data grid Z containing elevations in units of
%   meters.  The two points are provided as vectors of latitudes and
%   longitudes in units of degrees.  The resulting logical variable VIS is
%   equal to one when the two points are visible to each other, and zero
%   when the line of sight is obscured by terrain.  If any of the  input
%   arguments are empty, LOS2 attempts to gather the data from the current
%   axes.  With one or more output arguments, no figures are created and
%   only the data is returned.  R can be a spatialref.GeoRasterReference
%   object, a referencing vector, or a referencing matrix.
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
%   VIS = LOS2(Z, R, LAT1, LON1, LAT2, LON2, ALT1) places the first
%   point at the specified altitude in meters above the surface (on a
%   tower, for instance).  This is equivalent to putting the point on a
%   tower.  If omitted, point 1 is assumed to be on the surface.  ALT1 may
%   be either a scalar or a vector with the same length as LAT1, LON1,
%   LAT2, and LON2.
%
%   VIS = LOS2(Z, R, LAT1, LON1, LAT2, LON2, ALT1, ALT2) places both
%   points at a specified altitudes in meters above the surface.  ALT2 may
%   be either a scalar or a vector with the same length as LAT1, LON1,
%   LAT2, and LON2.  If ALT2 is omitted, point 2 is assumed to be on the
%   surface.
%
%   VIS = LOS2(Z, R, LAT1, LON1, LAT2, LON2, ALT1, ALT2, ALT1OPT)
%   controls the interpretation of ALT1 as either a relative altitude
%   (ALT1OPT equals 'AGL', the default) or an absolute altitude (ALT1OPT
%   equals 'MSL').  If the altitude option is 'AGL', ALT1 is interpreted as
%   the altitude of point 1 in meters above the terrain ("above ground
%   level").  If ALT1OPT is 'MSL', ALT1 is interpreted as altitude above
%   zero ("mean sea level").
%
%   VIS = LOS2(Z, R, LAT1, LON1, LAT2, LON2, ALT1, ALT2, ALT1OPT, ALT2OPT)
%   controls the interpretation ALT2. 
%
%   VIS = LOS2(Z, R, LAT1, LON1, LAT2, LON2, ALT1, ALT2, ALT1OPT,...
%   ALT2OPT, ACTUALRADIUS) does the visibility calculation on a sphere with
%   the specified radius.  If omitted, the radius of the earth in meters is
%   assumed.  The altitudes, elevations and the radius should be in the
%   same units.  This calling form is most useful for computations on
%   bodies other than the earth.
%
%   VIS = LOS2(Z, R, LAT1, LON1, LAT2, LON2, ALT1, ALT2, ALT1OPT,...
%   ALT2OPT, ACTUALRADIUS, EFFECTIVERADIUS) assumes a larger radius for
%   propagation of the line of sight.  This can account for the curvature
%   of the signal path due to refraction in the atmosphere.  For example,
%   radio propagation in the atmosphere is commonly treated as straight
%   line propagation on a sphere with 4/3rds the radius of the earth.  In
%   that case the last two arguments would be R_e and 4/3*R_e, where R_e is
%   the radius of the earth.  Use Inf as the effective radius for flat  
%   earth visibility calculations.  The altitudes, elevations and the radii
%   should be in the same units. 
%
%   [VIS, VISPROFILE, DIST ,H, LATTRK, LONTRK] = LOS2(...), for scalar
%   inputs (LAT1, LON1, etc.), returns vectors of points along the path
%   between the two points.  VISPROFILE is a logical vector containing true
%   (logical(1) where the intermediate points are visible and false
%   (logical(0)) otherwise.  DIST is the distance along the path (in meters
%   or the units of the radius).  H contains the terrain profile relative
%   to the vertical datum along the path.  LATTRK and LONTRK are the
%   latitudes and longitudes of the points along the path.  For vector
%   inputs LOS2 returns VISPROFILE, DIST, H, LATTRK, and LONTRK as cell
%   arrays, with one cell per element of LAT1, LON1, etc.
%
%   LOS2(...), with no output arguments, displays the visibility profile
%   between the two points in a new figure. 
%
%   Example
%   -------
%   Z = 500*peaks(100);
%   refvec = [1000 0 0];
%   [lat1, lon1, lat2, lon2] = deal(-0.027, 0.05, -0.093, 0.042);
%   los2(Z,refvec,lat1,lon1,lat2,lon2,100);
% 
%   figure;
%   axesm('globe','geoid',earthRadius('meters'))
%   meshm(Z, refvec, size(Z), Z); axis tight
%   camposm(-10,-10,1e6); camupm(0,0)
%   demcmap('inc', Z, 1000); shading interp; camlight
% 
%   [vis,visprofile,dist,h,lattrk,lontrk] = los2(Z,refvec,lat1,lon1,lat2,lon2,100);
%   plot3m(lattrk([1;end]),lontrk([1; end]),h([1; end])+[100; 0],'r','linewidth',2)
%   plotm(lattrk(~visprofile),lontrk(~visprofile),h(~visprofile),'r.','markersize',10)
%   plotm(lattrk(visprofile),lontrk(visprofile),h(visprofile),'g.','markersize',10)
% 
%   See also VIEWSHED.

% Copyright 1996-2010 The MathWorks, Inc.
% $Revision: 1.5.4.16 $  $Date: 2010/11/17 11:24:20 $
% Written by Walter Stumpf

if nargin < 2 || (isempty(map) && isempty(R))
   [map, R] = getrmm;
end

if nargin < 4 || (isempty(lat1) && isempty(lon1))
   disp('Click on the map for point 1')
   [lat1,lon1] = inputm(1);
end

if nargin < 6 || (isempty(lat2) && isempty(lon2))
   disp('Click on the map for point 2')
   [lat2,lon2] = inputm(1);
end

if nargin < 7; oalt = 100*eps; end % observer on the surface
if nargin < 8; talt = 0; end % look at terrain, not above it
if nargin < 9; oaltopt = 'AGL'; end % observer altitude above ground level
if nargin < 10; taltopt = 'AGL'; end % target above ground level
if nargin < 11; actualradius = earthRadius; end
if nargin < 12; apparentradius = actualradius; end % use Inf for flat earth LOS calculations

checklatlon(lat1, lon1, mfilename, 'LAT1', 'LON1', 3, 4)
checklatlon(lat2, lon2, mfilename, 'LAT2', 'LON2', 5, 6)

if numel(oalt) == 1
    oalt = repmat(oalt,size(lat1));
end

if numel(talt) == 1
    talt = repmat(talt,size(lat1));
end

% If R is already spatial referencing object, validate it. Otherwise
% convert the input referencing vector or matrix.
R = internal.map.convertToGeoRasterRef( ...
    R, size(map), 'degrees', mfilename, 'R', 2);

% loop over pairs of observer and target locations
if numel(lat1) == 1
    % Single pair of points:
    %    Output 1 is a logical scalar
    %    Output 2 is a logical array
    %    Outputs 3-6 are numerical arrays
    
    % compute the elevation profile between the start and end point
    [z, dist, lattrk, lontrk] = los2MapProfile( ...
        map, R, [lat1; lat2], [lon1; lon2], actualradius);

    % compute visibility of points on the track between the start and end point
    visprofile = losprofile(dist, z, oalt, talt, ...
                               oaltopt, taltopt, apparentradius, nargout);

    visprofile = reshape(visprofile,size(lattrk));
    
    vis = visprofile(end);
else
    % Multiple pairs of points:
    %   Output 1 is a logical array
    %   Outputs 2-6 are cell arrays that contain
    %      logical or numerical arrays
    vis = false(1, numel(lat1));
    visprofile = cell(1, numel(lat1));
    z      = cell(1, numel(lat1));
    dist   = cell(1, numel(lat1));
    lattrk = cell(1, numel(lat1));
    lontrk = cell(1, numel(lat1));
    
    for i = 1:length(lat1)
        % compute the elevation profile between the start and end point
        [z{i}, dist{i} ,lattrk{i}, lontrk{i}] = los2MapProfile(...
            map, R, [lat1(i); lat2(i)], [lon1(i); lon2(i)], actualradius);

        % compute visibility of points on the track between the start and end point
        visprofile{i} = losprofile(dist{i}, z{i}, oalt(i), talt(i), ...
            oaltopt, taltopt, apparentradius, nargout);

        visprofile{i} = reshape(visprofile{i},size(lattrk{i}));

        vis(i) = visprofile{i}(end);
    end
end

%--------------------------------------------------------------------------

function [z, rng, lat, lon] = los2MapProfile(map, R, lat, lon, radius)
% Core computations performed by MAPPROFILE and LTLN2VAL,
% streamlined for use within LOS2.

% cell extent in latitude
dlat = abs(R.DeltaLatNumerator / R.DeltaLatDenominator);

% interpolate points to less than the elevation grid spacing
[lat,lon] = doInterpm(lat, lon, 0.9*dlat, 'gc', 'degrees');

% The following 8 lines are equivalent to this:
%
%   z = interpGeoRaster(map, R, lat, lon, 'linear');
%
% but with a little less overhead.

% Interpolate bilinearly in intrinsic coordinates.
xi = R.longitudeToIntrinsicX(lon);
yi = R.latitudeToIntrinsicY(lat);

% Snap in all points that fall within distance 0.5 of an edge, so that
% we get a non-NaN value for them from interp2.
xi(0.5 <= xi & xi < 1) = 1;
yi(0.5 <= yi & yi < 1) = 1;

[M,N] = size(map);
xi(N < xi & xi <= N + 0.5) = N;
yi(M < yi & yi <= M + 0.5) = M;

z = interp2(map, xi, yi, '*bilinear');

dist = greatcircleinv( lat(1:end-1)*pi/180, lon(1:end-1)*pi/180, ...
    lat(2:end)*pi/180, lon(2:end)*pi/180, radius(1));

rng = [0; cumsum(dist)];

%--------------------------------------------------------------------------

function vis = losprofile(...
    rng, zin, oalt, talt, oaltopt, taltopt, apparentradius, nout)

rng = rng(:)';
zin = zin(:)';

if ~isinf(apparentradius)
    [x, z] = adjustterrain(rng, zin, apparentradius);
else
    x = rng;
    z = zin;
end

% Convert AGL observer altitude to MSL 
oaltopt = validatestring(oaltopt, {'AGL','MSL'}, 'LOS2', 'ALT1OPT', 9);
if strcmp(oaltopt,'AGL')
    %  Observer is at first location
    oalt =  z(1) + oalt;
end

% Shift terrain so observer is at altitude 0, and terrain altitudes are relative
% to the observer

z = z - oalt;  % Observer is at first location

% Compute the angles of sight from the observer to each point on the profile.
% measured positive up from the center of the sphere

ang = pi + atan2(z,x);
if x(1) == 0 && z(1) == 0
    ang(1) = pi/2;  % Look straight down at observer's location
end

% Find the cumulative maximum:  maxtohere(k) equals max(ang(1:k))
maxangtohere = cummax(ang);

% Adjust the angles for the altitude of the target height above ground level 
% or sea level and redo calculation of angles. This makes the obscuring factor
% the terrain only, and not any target height. To model stuff above the terrain 
% like a forest canopy, pass in a z vector that has the added heights.

taltopt = validatestring(taltopt, {'AGL','MSL'}, 'LOS2', 'ALT2OPT', 10);
switch taltopt
    case 'AGL'        
        if ~isinf(apparentradius)
            [x2, z2] = adjustterrain(rng, zin + talt, apparentradius);
            z2 = z2 - oalt;
        else
            z2 = z + talt;
            x2 = x;
        end        
    case 'MSL'
        if ~isinf(apparentradius)
            [x2, z2] = adjustterrain(rng, talt + zeros(size(zin)), apparentradius);
            z2 = z2 - oalt;
        else
            z2 = (talt)*ones(size(zin)) - oalt;
            x2 = x;
        end
end

% Compute line of sight angles again.

ang2 = pi + atan2(z2,x2);
if x2(1) == 0 && z2(1) == 0
    ang2(1) = pi/2;  % Look straight down at observer's location
end

% Visible are points that rise above the maximum angle of LOS of intervening 
% terrain.

vis = (ang2 >= maxangtohere);

% Visibility of first point below terrain need a special test, since
% it always passes the "angles of terrain up to here" test

if (z2(1) < z(1)) && (z(1) < 0)
    vis(1) = false;
end

% Display calculation if no output arguments in main function

if nout == 0
    h = NaN(4,1);

    figure;
    hold on
    h(1) = plot(x,z,'k');

    if vis(end)
        h(5) = plot([0;x2(end)],[0;z2(end)],'g');
    else
        h(5) = plot([0;x2(end)],[0;z2(end)],'r');
    end

    if any(vis)
        h(2) = plot(x2(vis),z2(vis),'g+');
    end

    if any(~vis)
        h(3) = plot(x2(~vis),z2(~vis),'ro');
    end

    h(4) = plot(0,0,'mx');
    axis equal

    labels = {'Terrain','Visible','Obscured','Observer','Line of Sight'};
    indx = find(~isnan(h));
    legend(h(indx),labels(indx))
    xlabel 'Horizontal Distance from Observer'
    ylabel 'Vertical Distance from Observer'
end

%-----------------------------------------------------------------------

function [x, z] = adjustterrain(rng, zin, apparentradius)

% Adjust the terrain slice for the curvature of the sphere. The radius
% may potentially be different from the actual body, for example to
% model refraction of radio waves.

% [z, y, x] = geodetic2ecef(rng/apparentradius, rng*0, zin, [apparentradius 0]);
% z = z - apparentradius;

r = apparentradius + zin;
phi = rng/apparentradius;
x = r .* sin(phi);
z = r .* cos(phi) - apparentradius;

%--------------------------------------------------------------------------

function y = cummax(x)
% Find the cumulative maximum of row vector x:
%
%     y(k) equals max(x(1:k)) for k in [1 N]
%
% where N is numel(x).

% Construct a matrix whose k-th row contains x(1,1:k) followed
% by N-k zeros.
A = tril(x(ones(numel(x),1),:));

% Find the largest value in each column of A; return result in a row vector.
y = max(A,[],2)';
