function h = worldmap(varargin)
%WORLDMAP Construct map axes for given region of world
%
%   WORLDMAP REGION or WORLDMAP(REGION) sets up an empty map axes with
%   projection and limits suitable to the part of the world specified in
%   REGION.  REGION may be a string or a cell array of strings. Permissible
%   strings include names of continents, countries, and islands, as well as
%   'World', 'North Pole', 'South Pole', and 'Pacific'.
%
%   WORLDMAP with no arguments presents a menu from which you can select
%   the name of a single continent, country, island, or other region.
%   
%   WORLDMAP(LATLIM, LONLIM) allows you to define a custom geographic
%   region in terms of its latitude and longitude limits in degrees. LATLIM
%   and LONLIM are two-element vectors of the form [southern_limit
%   northern_limit] and [western_limit eastern_limit], respectively.
%   
%   WORLDMAP(Z, R) derives the map limits from the extent of a regular data
%   grid georeferenced by R.  R can be a spatialref.GeoRasterReference
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
%   H = WORLDMAP(...) returns the handle of the map axes.
%  
%   For cylindrical projections, WORLDMAP uses TIGHTMAP set the axis limits
%   tight around the map. If you change the projection, or just want more
%   white space around the map frame, use TIGHTMAP again or AXIS AUTO.
%
%   Example 1
%   ---------
%   % World map with coarse coastlines
%   worldmap('World')
%   load coast
%   plotm(lat, long)
%
%   Example 2
%   ---------
%   % Worldmap with land areas, major lakes and rivers, and cities and
%   % populated places
%   ax = worldmap('World');
%   setm(ax, 'Origin', [0 180 0])
%   land = shaperead('landareas', 'UseGeoCoords', true);
%   geoshow(ax, land, 'FaceColor', [0.5 0.7 0.5])
%   lakes = shaperead('worldlakes', 'UseGeoCoords', true);
%   geoshow(lakes, 'FaceColor', 'blue')
%   rivers = shaperead('worldrivers', 'UseGeoCoords', true);
%   geoshow(rivers, 'Color', 'blue')
%   cities = shaperead('worldcities', 'UseGeoCoords', true);
%   geoshow(cities, 'Marker', '.', 'Color', 'red')
%
%   Example 3
%   ---------
%   % Map of Antarctica
%   worldmap('antarctica')
%   antarctica = shaperead('landareas', 'UseGeoCoords', true,...
%       'Selector',{@(name) strcmp(name,'Antarctica'), 'Name'});
%   patchm(antarctica.Lat, antarctica.Lon, [0.5 1 0.5])
%
%   Example 4
%   ---------
%   % Map of Africa and India with major cities and populated places
%   worldmap({'Africa','India'})
%   land = shaperead('landareas.shp', 'UseGeoCoords', true);
%   geoshow(land, 'FaceColor', [0.15 0.5 0.15])
%   cities = shaperead('worldcities', 'UseGeoCoords', true);
%   geoshow(cities, 'Marker', '.', 'Color', 'blue')
%
%   Example 5
%   ---------
%   % Map of the geoid over South America and the central Pacific
%   worldmap([-50 50],[160 -30])
%   load geoid
%   geoshow(geoid, geoidrefvec, 'DisplayType', 'texturemap');
%   load coast
%   geoshow(lat, long)
%
%   Example 6
%   ---------
%   % Map of terrain elevations in Korea
%   load korea
%   h = worldmap(map, refvec);
%   set(h, 'Visible', 'off')
%   geoshow(h, map, refvec, 'DisplayType', 'texturemap')
%   demcmap(map)
%
%   Example 7
%   ---------
%   % Map of the United States of America
%   ax = worldmap('USA');
%   load coast
%   geoshow(ax, lat, long,...
%       'DisplayType', 'polygon', 'FaceColor', [.45 .60 .30])
%   states = shaperead('usastatelo', 'UseGeoCoords', true);
%   faceColors = makesymbolspec('Polygon',...
%       {'INDEX', [1 numel(states)], 'FaceColor', polcmap(numel(states))});
%   geoshow(ax, states, 'DisplayType', 'polygon', 'SymbolSpec', faceColors)
%
%   See also AXESM, FRAMEM, GEOSHOW, GRIDM, MLABEL, PLABEL, TIGHTMAP, USAMAP.

% Copyright 1996-2011 The MathWorks, Inc.
% $Revision: 1.11.4.6.2.1 $  $Date: 2011/01/29 14:47:38 $

ax = regionmap(mfilename, varargin{:});

% Avoid command-line output if no output variable is specified.
if nargout == 1
    h = ax;
end
