function [vismap, R] = viewshed(map, R, lat1, lon1, ...
    oalt, talt, oaltopt, taltopt, actualradius, apparentradius)
%VIEWSHED Areas visible from point on terrain elevation grid
%
%   [VIS, R] = VIEWSHED(Z, R, LAT, LON) computes areas visible from a
%   point on a digital elevation grid, Z. The elevations are provided as
%   a regular data grid containing elevations in units of meters. The
%   observer location is provided as scalar latitude and longitude in
%   units of degrees. The visibility grid VIS contains ones at the
%   surface locations visible from the observer location, and zeros
%   where the line of sight is obscured by terrain.  R can be a
%   spatialref.GeoRasterReference object, a referencing vector, or a
%   referencing matrix.
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
%                   [lon lat] = [row col 1] * R.
%
%   If R is a referencing matrix, it must define a (non-rotational,
%   non-skewed) relationship in which each column of the data grid falls
%   along a meridian and each row falls along a parallel. The value of R
%   on output is identical to the value supplied as input. The observer
%   should be located within the latitude-longitude limits of the
%   elevation grid. If the observer is located outside the grid, there
%   is insufficient information to calculate a viewshed.  In this case
%   VIEWSHED issues a warning and sets all elements of VIS to zero.
% 
%   VIEWSHED(Z, R, LAT, LON, observerAltitude) places the observer at
%   the specified altitude in meters above the surface. This is
%   equivalent to putting the observer on a tower. If omitted, the
%   observer is assumed to be on the surface.
% 
%   VIEWSHED(Z, R, LAT, LON, observerAltitude, targetAltitude) checks
%   for visibility of target points a specified distance above the
%   terrain. This is equivalent to putting the target points on towers,
%   but the towers do not obstruct the view. if omitted, the target
%   points are assumed to be on the surface.
% 
%   VIEWSHED(Z, R, LAT, LON, observerAltitude, targetAltitude, ...
%     observerAltitudeOption) controls whether the observer is at a
%   relative or absolute altitude. If observerAltitudeOption is 'AGL',
%   then observerAltitude is in meters above ground level. If
%   observerAltitudeOption is 'MSL', observerAltitude is interpreted as
%   altitude above zero, or mean sea level. If omitted, 'AGL' is
%   assumed.
% 
%   VIEWSHED(Z, R, LAT, LON, observerAltitude, targetAltitude, ...
%      observerAltitudeOption, targetAltitudeOptionOption) controls
%   whether the target points are at a relative or absolute altitude. If
%   the target altitude option is 'AGL', then targetAltitude is in
%   meters above ground level. If targetAltitudeOption is 'MSL', then
%   targetAltitude is interpreted as altitude above zero, or mean sea
%   level. If omitted, 'AGL' is assumed.
% 
%   VIEWSHED(Z, R, LAT, LON, observerAltitude, targetAltitude, ...
%      observerAltitudeOption, targetAltitudeOption, actualRadius)
%   does the visibility calculation on a sphere with the specified
%   radius. If omitted, the radius of the earth in meters is assumed.
%   The altitudes, elevations and the radius should be in the same
%   units. This calling form is most useful for computations on bodies
%   other than the earth.
% 
%   VIEWSHED(Z, R, LAT, LON, observerAltitude, targetAltitude, ...
%      observerAltitudeOption, targetAltitudeOption, actualRadius, ...
%      effectiveRadius) assumes a larger radius for propagation of the
%   line of sight. This can account for the curvature of the signal path
%   due to refraction in the atmosphere. For example, radio propagation
%   in the atmosphere is commonly treated as straight line propagation
%   on a sphere with 4/3rds the radius of the earth. In that case the
%   last two arguments would be R_e and 4/3*R_e, where R_e is the radius
%   of the earth. Use Inf for flat earth viewshed calculations. The
%   altitudes, elevations and the radii should be in the same units. 
% 
%   Example
%   -------
%    Z = 500*peaks(100);
%    refvec = [ 1000 0 0];
%    [lat1,lon1,lat2,lon2]=deal(-0.027,0.05,-0.093,0.042);
%    
%    [visgrid,visleg] = viewshed(Z,refvec,lat1,lon1,100);
%    [vis,visprofile,dist,zi,lattrk,lontrk] ...
%        = los2(Z,refvec,lat1,lon1,lat2,lon2,100);
%  
%    axesm('globe','geoid',earthRadius('meters'))
%    meshm(visgrid,visleg,size(Z),Z); axis tight
%    camposm(-10,-10,1e6); camupm(0,0)
%    colormap(flipud(summer(2))); brighten(0.75);
%    shading interp; camlight
%    h = lcolorbar({'obscured','visible'});
%    set(h,'Position',[.875 .45 .02 .1])
%  
%    plot3m(lattrk([1;end]),lontrk([1; end]),zi([1; end])+[100; 0],'r','linewidth',2)
%    plotm(lattrk(~visprofile),lontrk(~visprofile),zi(~visprofile),'r.','markersize',10)
%    plotm(lattrk(visprofile),lontrk(visprofile),zi(visprofile),'g.','markersize',10)
%
%   See also LOS2.

% Copyright 1996-2010 The MathWorks, Inc.
% $Revision: 1.5.4.15 $  $Date: 2010/11/17 11:24:39 $
% Written by Walter Stumpf and Betty Youmans.

% Elevations assumed to be in meters on the earth. If not on the earth,
% provide the actual (and optionally apparent) radius of the spherical
% body in the same units as the elevation data.

error(nargchk(4, Inf, nargin, 'struct'))

if nargin < 5; oalt = 10*eps; end % observer on the surface
if nargin < 6; talt = 0; end % look at terrain, not above it
if nargin < 7; oaltopt = 'AGL'; end % observer altitude above ground level
if nargin < 8; taltopt = 'AGL'; end % target above ground level
if nargin < 9; actualradius = earthRadius; end 
if nargin < 10; apparentradius = actualradius; end; % use Inf for flat earth LOS calculations

checklatlon(lat1, lon1, mfilename, 'LAT', 'LON', 3, 4)

% If R is already spatial referencing object, validate it. Otherwise
% convert the input referencing vector or matrix. (Save the original
% value to use as the second output argument.)
R_input = R;
R = internal.map.convertToGeoRasterRef( ...
    R, size(map), 'degrees', 'VIEWSHED', 'R', 2);

oaltopt = upper(oaltopt);
taltopt = upper(taltopt);

[nr,nc] = size(map);
vismap = NaN(nr,nc);

% Throw warning if observer is outside the grid.
if ~ingeoquad(lat1, lon1, R.Latlim, R.Lonlim)
    warning('map:viewshed:observerOutsideGrid', ...
        ['The observer is located outside the terrain grid;\n' ...
         'visibility will be set to 0 (false) for all grid cells'])
end

for i=1:nr
   [lat2, lon2] = R.intrinsicToGeographic(1,i);
   
   [~,vis,~,~,lat,lon] = los2(map,R,lat1,lon1,lat2,lon2,...
      oalt,talt,oaltopt,taltopt,actualradius,apparentradius);
     
   vismap = embed(lat,lon,vis,vismap,R);
      
   [lat2,lon2] = R.intrinsicToGeographic(nc,i);
      
   [~,vis,~,~,lat,lon] = los2(map,R,lat1,lon1,lat2,lon2,...
      oalt,talt,oaltopt,taltopt,actualradius,apparentradius);
   
   vismap = embed(lat,lon,vis,vismap,R);   
end


for i=1:nc
   [lat2,lon2] = R.intrinsicToGeographic(i,1);
  
   [~,vis,~,~,lat,lon] = los2(map,R,lat1,lon1,lat2,lon2,...
      oalt,talt,oaltopt,taltopt,actualradius,apparentradius);
   
   vismap = embed(lat,lon,vis,vismap,R);
   
   [lat2,lon2] = R.intrinsicToGeographic(i,nr);
   
   [~,vis,~,~,lat,lon] = los2(map,R,lat1,lon1,lat2,lon2,...
      oalt,talt,oaltopt,taltopt,actualradius,apparentradius);
   
   vismap = embed(lat,lon,vis,vismap,R);
end

% Use los2 to handle any remaining NaNs in vismap

nanvec=isnan(vismap);
[inan,jnan]=find(nanvec);

for i=1:length(inan),
    [lat2,lon2] = R.intrinsicToGeographic(jnan(i),inan(i));
    
    vis = los2(map,R,lat1,lon1,lat2,lon2,...
        oalt,talt,oaltopt,taltopt,actualradius,apparentradius);
    
    vismap = embed(lat2,lon2,vis,vismap,R);
end

R = R_input;

%--------------------------------------------------------------------------

function Z = embed(lat, lon, value, Z, R)
%EMBED  Encode data points into regular data grid
%
%   Streamlined version of the public IMBEDM function
%
%   Z = EMBED(LAT, LON, VALUE, Z, R) resets certain entries of a regular
%   data grid, Z.  R is a spatialref.GeoRasterReference object.  Its
%   RasterSize property must be consistent with size(Z). The entries to be
%   reset correspond to the locations defined by the latitude and longitude
%   position vectors LAT and LON. The entries are reset to the same number
%   if VALUE is a scalar, or to individually specified numbers if VALUE is
%   a vector the same size as LAT and LON. If any points lie outside the
%   input grid, a warning is issued.  All input angles are in degrees.

%  Eliminate NaNs from the input data
qNaN = isnan(lat) | isnan(lon);
lat(qNaN) = [];
lon(qNaN) = [];
value(qNaN) = [];

%  Convert the lat and lon data to cell positions
[r, c] = R.geographicToSub(lat, lon);

% Remove any points that fall outside the grid
qOutside = isnan(r);
r(qOutside) = [];
c(qOutside) = [];
value(qOutside) = [];

%  Embed the points into the grid
indx = (c-1)*size(Z,1) + r;
Z(indx) = value;
