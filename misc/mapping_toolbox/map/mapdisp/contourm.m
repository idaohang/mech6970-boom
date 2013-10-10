function [c, g] = contourm(varargin)
%CONTOURM Project 2-D contour plot of map data
%
%  CONTOURM(Z,R) creates a contour plot of the regular data grid Z.
%  R can be a spatialref.GeoRasterReference object, a referencing
%  vector, or a referencing matrix.
%
%  If R is a spatialref.GeoRasterReference object, its RasterSize
%  property must be consistent with size(Z).
%
%  If R is a referencing vector, it must be a 1-by-3 with elements:
%
%     [cells/degree northern_latitude_limit western_longitude_limit]
%
%  If R is a referencing matrix, it must be 3-by-2 and transform raster
%  row and column indices to/from geographic coordinates according to:
%
%                    [lon lat] = [row col 1] * R.
%
%  If R is a referencing matrix, it must define a (non-rotational,
%  non-skewed) relationship in which each column of the data grid falls
%  along a meridian and each row falls along a parallel. If the current
%  axis is a map axis, the coordinates of Z will be projected using the
%  projection structure from the axis. The contours are drawn at their
%  corresponding Z level.
%
%  CONTOURM(LAT, LON, Z) displays a contour plot of the geolocated M-by-N
%  data grid, Z.  LAT and LON can be the size of Z or can specify the 
%  corresponding row and column dimensions for Z.
%
%  CONTOURM(Z, R, N) or CONTOURM(LAT, LON, Z, N) where N is a positive
%  scalar integer, draws N contour levels.
%
%  CONTOURM(Z, R, V) or CONTOURM(LAT, LON, Z, V) where V is a vector, draws 
%  contours at the levels specified by the input vector V. Use V = [v v] to 
%  compute a single contour at level v. 
%
%  CONTOURM(..., LINESPEC) uses any valid LineSpec string to draw the 
%  contour lines.
%
%  CONTOURM(..., PARAM1, VAL1, PARAM2, VAL2, ...) allows you to set the
%  following optional parameters: Fill, LevelStep, LineColor, LineStyle,
%  LineWidth, and ShowText. See the CONTOURM reference page for full
%  descriptions. In addition, any of the following hggroup properties
%  may be specified: HandleVisibility, Parent, Tag, UserData, and
%  Visible.
%
%  C = CONTOURM(...) returns a standard contour matrix, C, with the first
%  row representing longitude data and the second row representing latitude
%  data.
%
%  [C,H] = CONTOURM(...) returns the contour matrix and the handle to the
%  contour patches drawn onto the current axes. The handle is type hggroup.
%
%  % Example 1
%  % ---------
%  % Contour the EGM96 geoid heights, label them, and add a legend.
%  load geoid
%  figure
%  [c,h] = contourm(geoid,geoidrefvec,'LevelStep',20,'ShowText','on');
%  xlabel('Longitude')
%  ylabel('Latitude')
%  clegendm(c,h,-1)
%
%  % Example 2
%  % ---------
%  % Contour geoid heights for an area including Korea with a backdrop of
%  % terrain elevations and bathymetry.
%
%  % Load the data.
%  load korea
%  load geoid
%
%  % Create a map axes that includes Korea.
%  figure
%  worldmap(map, refvec)
%
%  % Display the digital elevation data and apply a colormap.
%  geoshow(map, refvec, 'DisplayType', 'texturemap');
%  demcmap(map)
%
%  % Contour the geoid values from -100 to 100 in increments of 5.
%  [c,h] = contourm(geoid, geoidlegend, -100:5:100, 'k');
%
%  % Add red labels with white backgrounds to the contours.
%  ht = clabelm(c,h);
%  set(ht,'Color','r','BackgroundColor','white','FontWeight','bold')
%
%  See also CLABELM, CONTOUR, CONTOUR3M, CONTOURFM, GEOSHOW.

% Copyright 1996-2010 The MathWorks, Inc.
% $Revision: 1.3.4.13 $  $Date: 2010/11/17 11:25:07 $

% Catch deprecated syntax.
if nargin == 0
    errorOnMissingUI('contourm')
end

error(nargchk(2, inf, nargin, 'struct'))

% Get the data grid, reference object, contour level list vector,
% and a cell array of optional property name-value pairs.
[Z, R, levelList, pvpairs] = parseContourInputs(varargin);

% Obtain axes, and get ready to plot if it's a map axes.
[ax, pvpairs] = internal.map.findNameValuePair('Parent', gca, pvpairs{:});
if ismap(ax)
    nextmap(ax)
end

% Construct a GeoContourGroup object and associated hggroup
h = internal.mapgraph.GeoContourGroup(ax, Z, R, levelList);

% Set any properties supplied by the user (or contourfm or contour3m).
if ~isempty(pvpairs)
    set(h, pvpairs{:})
end

% Construct the contour lines and polygons themselves.
h.refresh()

% Set the Tag property of the hggroup unless it's already been set
% (which could happen if it appeared in the property-value list).
tag = get(h.HGGroup,'Tag');
if isempty(tag)
    set(h.HGGroup,'Tag','Contour')
end

if nargout > 0
    c = geostructToContourMatrix(h.getContourLines());
end

if nargout > 1
    g = h.HGGroup;
end

%-------------------------------------------------------------------

function c = geostructToContourMatrix(L)
% Convert contour line geostruct L to geographic contour matrix c.

% Allocate contour matrix.
ncols = 0;
for k = 1:numel(L)
    ncols = ncols + numel(L(k).Lon) + 1;
end
c = zeros(2,ncols);

% Fill in contour matrix.
n = 1;
for k = 1:numel(L)
    % Process the k-th level.
    lonk = L(k).Lon;
    latk = L(k).Lat;
    [first, last] = internal.map.findFirstLastNonNan(lonk);
    for j = 1:numel(first)
        % Process the j-th part of the k-th level.
        s = first(j);
        e = last(j);
        lon = lonk(s:e);
        lat = latk(s:e);
        count = numel(lon);
        c(:,n) = [L(k).Level; count];
        m = n + count;
        n = n + 1;
        c(1,n:m) = lon(:)';
        c(2,n:m) = lat(:)';
        n = m + 1;
    end
end
c(:,n:end) = [];
