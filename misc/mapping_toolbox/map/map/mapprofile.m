function [zout,rng,lat,lon] = mapprofile(varargin)
%MAPPROFILE Interpolate between waypoints on regular data grid
%
% MAPPROFILE plots a profile of values between waypoints on 
% a displayed regular data grid. MAPPROFILE uses the current
% object if it is a regular data grid, or the first regular 
% data grid found on the current axes. The grid's zdata is 
% used for the profile. The color data is used in the absence 
% of zdata. The result is displayed in a new figure.
%
% [zi,rng,lat,lon] = MAPPROFILE returns the values of the profile
% without displaying them. The output zi contains interpolated 
% values from map along great circles between the waypoints. 
% rng is a vector of associated distances from the first waypoint 
% in units of degrees of arc along the surface. lat and lon are 
% the corresponding latitudes and longitudes. 
%
% [zi,rng,lat,lon] = MAPPROFILE(Z,R,lat,lon) accepts as input a regular
% data grid and waypoint vectors. No displayed grid is required. Sets of
% waypoints may be separated by NaNs into line sequences. The output
% ranges are measured from the first waypoint within a sequence.  R can
% be a spatialref.GeoRasterReference object, a referencing vector, or a
% referencing matrix.
%
% If R is a spatialref.GeoRasterReference object, its RasterSize
% property must be consistent with size(Z).
%
% If R is a referencing vector, it must be a 1-by-3 with elements:
%
%     [cells/degree northern_latitude_limit western_longitude_limit]
%
% If R is a referencing matrix, it must be 3-by-2 and transform raster
% row and column indices to/from geographic coordinates according to:
% 
%                  [lon lat] = [row col 1] * R.
%
% If R is a referencing matrix, it must define a (non-rotational,
% non-skewed) relationship in which each column of the data grid falls
% along a meridian and each row falls along a parallel.
%
% [zi,rng,lat,lon] = MAPPROFILE(Z,R,lat,lon,rngunits)
% specifies the units of the output ranges along the profile. 
% Valid range units inputs are any distance string recognized by
% UNITSRATIO. Surface distances are computed using the default 
% radius of the earth. If omitted, 'degrees' are assumed.
%
% [zi,rng,lat,lon] = MAPPROFILE(Z,R,lat,lon,ellipsoid) 
% uses the provided ellipsoid definition in computing the range
% along the profile. The ellipsoid vector is of the form
% [semimajor axes, eccentricity].  The output range is reported in
% the same distance units as the semimajor axes of the ellipsoid
% vector. If omitted, the range vector is for a sphere.
%
% [zi,rng,lat,lon] = MAPPROFILE(Z,R,lat,lon,rngunits,...
% 'trackmethod','interpmethod') and
% [zi,rng,lat,lon] = MAPPROFILE(Z,R,lat,lon,ellipsoid,...
% 'trackmethod','interpmethod') control the interpolation methods
% used. Valid trackmethods are 'gc' for great circle tracks 
% between waypoints, and 'rh' for rhumb lines. Valid interpmethods
% for interpolation within the matrix are 'bilinear' for linear 
% interpolation, 'bicubic' for cubic interpolation, and 'nearest' 
% for nearest neighbor interpolation. If omitted, 'gc' and 'bilinear'
% are assumed.
%
% See also LTLN2VAL, LOS2.

% Copyright 1996-2010 The MathWorks, Inc.
% $Revision: 1.5.4.12 $  $Date: 2010/11/17 11:24:23 $

% Defaults for optional arguments
rngunits = 'deg'; 
trackmethod = 'gc';
interpmethod = 'bilinear';

% other inputs
error(nargchk(0, 7, nargin, 'struct'))

if numel(varargin) == 0 % get from figure
   [map,R] = getrmm;
   [lat,lon] = inputm;
elseif length(varargin) == 4
   [map,R,lat,lon] = deal(varargin{1:4});
elseif length(varargin) == 5
   [map,R,lat,lon,rngunits] = deal(varargin{1:5});
elseif length(varargin) == 6
   [map,R,lat,lon,rngunits,trackmethod] = deal(varargin{1:6});
elseif length(varargin) == 7
   [map,R,lat,lon,rngunits,trackmethod,interpmethod] = deal(varargin{1:7});
end

% check if geoid was provided
if isnumeric(rngunits)
    geoid = rngunits;
    geoid = geoidtst(geoid);
else
    geoid = [];
end

checklatlon(lat, lon, mfilename, 'LAT', 'LON', 3, 4)

% check trackmethod
switch trackmethod
case {'gc','rh'}
   % OK
otherwise
   error(['map:' mfilename ':mapError'], ...
       'Recognized method strings are ''gc'' and ''rh''') 
end

%  Try to ensure vectors don't begin or end with NaNs
if isnan(lat(1)) || isnan(lon(1))
	lat = lat(2:end);
	lon = lon(2:end);
end
if isnan(lat(end)) || isnan(lon(end))
	lat = lat(1:end-1);
	lon = lon(1:end-1);
end

%  If R is already spatial referencing object, validate it. Otherwise
%  convert the input referencing vector or matrix.
R = internal.map.convertToGeoRasterRef( ...
    R, size(map), 'degrees', 'mapprofile', 'R', 2);

% determine distances along the track, starting from zero at the
% beginning of each segment.
[z, rng, lat, lon] = doMapProfile(...
    map, R, lat, lon, geoid, trackmethod, interpmethod);

% Convert ranges to desired distance units. If geoid provided,
% output is in units of geoid(1).
if isempty(geoid)
    rng = deg2dist(rng, rngunits);
end

%if no output arguments, plot results
if nargout == 0;   
   outputplot(nargin,numel(lat),lat,lon,rng,z,rngunits)
else
   zout = z;
end

%-----------------------------------------------------------------------

function outputplot(nin,npts,lat,lon,rng,z,rngunits)

% displays output in a new figure

if npts > 2 % plot on a map
   
   if nin == 0;
      
      % Display results on a partial copy of the original map (line data only)
      hax = gca;
      hline = handlem('allline');
      mstruct = getm(hax);
      figure
      axesm miller
      set(gca,'UserData',mstruct)
      copyobj(hline,gca);
      
   else
      
      % Display results on a new map
      latlim = [min(lat(:)) max(lat(:))];
      lonlim = [min(lon(:)) max(lon(:))];
      latlim = mean(latlim)+1.5*diff(latlim)*[-1 1];
      lonlim = mean(lonlim)+1.5*diff(lonlim)*[-1 1];
      figure
      worldmap(latlim,lonlim)
      
   end
   
   % plot additional elements on a map axes
   
   framem on; gridm on; mlabel on; plabel on
   
   zdatam('frame',0)
   zdatam('alltext',0)
   zdatam('allline',0)
   zdatam('grid',0)     
   
   plot3m(lat,lon,z)
   plotm(lat,lon,':')
   stem3m(lat(:),lon(:),z(:))
   
   tightmap
   box off
   
   view(3)
   set(gca,'DataAspectRatio', [ 1 1 5*diff(zlim)/max(diff(xlim),diff(ylim))])

else   
   
   % 2-d plot
   figure
   plot(rng,z)
   if ischar(rngunits);
      xlabel(['Range [' rngunits ']' ])
   end
   ylabel 'Value'
   
end

%--------------------------------------------------------------------------

function dist = deg2dist(deg,units)
% Convert from spherical distance in degrees
%
%   DIST = DEG2DIST(DEG, UNITS) converts distances from degrees, as
%   measured along a great circle on a sphere with a radius of 6371 km
%   (the mean radius of the Earth) to the UNITS of length or angle
%   specified by the string  in UNITS.  If UNITS is 'degrees' or
%   'radians', then DIST itself will be a spherical distance.

angleUnits = {'degrees','radians'};
k = find(strncmpi(deblank(units), angleUnits, numel(deblank(units))));
if numel(k) == 1
    % In case units is 'degrees' or 'radians'
    dist = fromDegrees(angleUnits{k}, deg);
else
    % Assume that units specifies a length unit; convert using
    % kilometers as an intermediate unit.
    dist = unitsratio(units,'km') * deg2km(deg);
end

%--------------------------------------------------------------------------

function [z, rng, lat, lon] = doMapProfile(...
    map, R, lat, lon, ellipsoid, trackmethod, interpmethod)
% Core computations performed by MAPPROFILE.
% R is a spatialref.GeoRasterReference object.

% interpolate points to less than the elevation grid spacing
dlat = abs(R.DeltaLatNumerator / R.DeltaLatDenominator);
[lat,lon] = doInterpm(lat, lon, 0.9*dlat, trackmethod, 'deg');

[latcells,loncells] = polysplit(lat,lon);

% extract the elevation profile

inGrid = ~isnan(lat);
z = NaN + zeros(size(lat));
z(inGrid) = interpGeoRaster(map, R, lat(inGrid), lon(inGrid), interpmethod);

rng = [];
for i=1:numel(latcells)
    leglat = latcells{i};
    leglon = loncells{i};
    
    if isempty(ellipsoid)
        legdist = distance(leglat(1:end-1),leglon(1:end-1),leglat(2:end),leglon(2:end));
    else % ellipsoid provided
        legdist = distance(leglat(1:end-1),leglon(1:end-1),leglat(2:end),leglon(2:end),ellipsoid);
    end
    
    legrng = [0; cumsum(legdist)];
    if i==1;
        rng = [rng; legrng];
    else
        rng = [rng; NaN;legrng];
    end
end
