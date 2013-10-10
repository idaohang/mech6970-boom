function [latb,lonb] = bufferm(lat,lon,bufwidth,direction,npts,outputformat)
%BUFFERM Buffer zones for latitude-longitude polygons
%
%   [LATB, LONB] = BUFFERM(LAT, LON, BUFWIDTH) computes the buffer zone
%   around a line or polygon. If the vectors LAT and LON, in units of
%   degrees, define a line, then LATB and LONB define a polygon that
%   contains all the points that fall within a certain distance,
%   BUFWIDTH, of the line. BUFWIDTH is a scalar specified in degrees of
%   arc along the surface. If the vectors LAT and LON define a polygon,
%   then LATB and LONB define a region that contains all the points
%   exterior to the polygon that fall within BUFWIDTH of the polygon.
%
%   [LATB, LONB] = BUFFERM(LAT, LON, BUFWIDTH, DIRECTION) uses the optional
%   string DIRECTION to specify whether the buffer zone is inside ('in')
%   or outside ('out') of the polygon. If you do not supply a direction
%   string, BUFFERM uses 'out' as the default and returns a buffer zone
%   outside the polygon. If you supply 'in' as the direction string,
%   BUFFERM returns a buffer zone inside the polygon. If you are finding
%   the buffer zone around a line, 'out' is the only valid option.
%
%   [LATB, LONB] = BUFFERM(LAT, LON, DIST, DIRECTION, NPTS) controls the
%   number of points used to construct circles about the vertices of the
%   polygon. A larger number of points produces smoother buffers, but
%   requires more time. If NPTS is omitted, 13 points per circle are
%   used.
%
%   Example
%   -------
%   load conus
%   tol = 0.05; % Tolerance for simplifying polygon outlines
%   [latr, lonr] = reducem(gtlakelat, gtlakelon, tol);
%   bufwidth = 1;  % Buffer width in degrees
%   [latb, lonb] = bufferm(latr, lonr, bufwidth, 'out');
%   [lati, loni] = bufferm(latr, lonr, 0.3*bufwidth, 'in');
%   figure('Color','w')
%   ax = usamap({'MN','NY'});
%   setm(ax,'MLabelLocation',5)
%   geoshow(latb, lonb, 'DisplayType', 'polygon', 'FaceColor', 'yellow')
%   geoshow(latr, lonr, 'DisplayType', 'polygon', 'FaceColor', 'blue')
%   geoshow(lati, loni, 'DisplayType', 'polygon', 'FaceColor', 'magenta')
%   geoshow(uslat, uslon)
%   geoshow(statelat, statelon)
%
%   See also POLYBOOL.

% Copyright 1996-2010 The MathWorks, Inc.
% $Revision: 1.5.4.13 $  $Date: 2010/06/26 04:57:08 $

% The following syntax remains supported but is not recommended and is
% intentionally undocumented:
%
%   [LATB, LONB] = BUFFERM(LAT, LON, BUFWIDTH, DIRECTION, NPTS, OUTPUTFORMAT)
%   controls the format of the returned buffer zones. OUTPUTFORMAT 'vector'
%   returns NaN-clipped vectors. OUTPUTFORMAT 'cutvector' returns
%   NaN-clipped vectors with cuts connecting holes to the exterior of the
%   polygon. OUTPUTFORMAT 'cell' returns cell arrays in which each element
%   of the cell array is a separate polygon. Each polygon may consist of
%   an outer contour followed by holes separated with NaNs.

% Validate argument BUFWIDTH.
validateattributes(bufwidth, {'numeric'}, ...
    {'positive','finite','scalar'}, mfilename, 'BUFWIDTH')

% Set/validate optional 4-th argument DIRECTION.
if nargin < 4
    direction = 'out';
else
    assert(any(strcmp(direction,{'in','out'})), ...
        'map:bufferm:unknownDirectionFlag', ...
        'Direction must be either ''%s'' or ''%s''.','in','out')
end

% Set/validate optional 5-th argument NPTS.
if nargin < 5
    npts = 13;
else
    validateattributes(npts, {'numeric'}, ...
        {'positive','integer','scalar'}, mfilename, 'NPTS')
    if mod(npts,2) == 0
        % Ensure that npts is odd
        npts = npts + 1;
    end
end

% Set/validate optional 6-th argument OUTPUTFORMAT.
if nargin < 6
    outputformat = 'vector';
else
    formats = {'vector','cell','cutvector'};
    assert(any(strcmp(outputformat,formats)), ...
        'map:bufferm:unknownFormatFlag', ...
        'Output format must be ''%s'', ''%s'', or ''%s''.', formats{:})        
end

% Validate/convert vertex arrays: convert possible cell array inputs to
% NaN-separated form; remove duplicate vertices; ensure column vectors.
if iscell(lat)
    checkcellvector(lat,lon)
    [lat, lon] = polyjoin(lat,lon);
else
    checklatlon(lat,lon,mfilename,'LAT','LON',1,2);
end
[lat, lon] = removeDuplicateVertices(lat, lon);
lat = lat(:);
lon = lon(:);

% Perform buffering operations using NaN-separated vectors.
[latb, lonb] = dobufferm(lat, lon, bufwidth, direction, npts - 1);

% Convert to alternative output format if requested.
if ~strcmp(outputformat,'vector')
    [latb, lonb] = polysplit(latb, lonb);
    if strcmp(outputformat,'cutvector')
        [latb,lonb] = polycut(latb,lonb);
    end
end

%-----------------------------------------------------------------------

function checkcellvector(lat, lon)

if ~isa(lat,'cell') || ~isa(lon,'cell')
    error('map:bufferm:expectedCellArrays', ...
        'Expected latitude and longitude to be cell arrays.')
end

if ~isequal(size(lat),size(lon))
    error('map:bufferm:inconsistentSizes', ...
        'Inconsistent dimensions on latitude and longitude input.');
end

for k = 1:numel(lat)
    if ~isequal(size(lat{k}),size(lon{k}))
        error('map:bufferm:inconsistentSizesInCell', ...
            'Inconsistent latitude and longitude dimensions within a cell.')
    end
    
    [lat{k}, lon{k}] = removeDuplicateVertices(lat{k}, lon{k});
end

%-----------------------------------------------------------------------

function [latb, lonb] = dobufferm(lat, lon, bufwidth, direction, n)
% Operate on NaN-separated column vectors with no duplicate vertices.

% Work with row vectors throughout, but keep track of input shape.
rowVectorInput = (size(lat,2) > 1);
lat = lat(:)';
lon = lon(:)';

latlim = [ -90  90];
lonlim = [-180 180];

% Identify open and closed curves.
[first, last] = internal.map.findFirstLastNonNan(lat);
isClosed = ((lat(first) == lat(last)) ...
    & (mod(lon(first) - lon(last), 360) == 0));
firstClosed = first(isClosed);
lastClosed  = last(isClosed);

if strcmp(direction,'in')
    % Process polygons only, as defined by closed curves.
    if ~isempty(firstClosed)
        % Buffer the edges of closed polygons using a planar topology;
        % convert the input polygons to the same topology and longitude
        % limits; keep only the buffered areas that fall within an
        % input polygon.
        [latClosed, lonClosed] = subsetCurves( ...
            lat, lon, firstClosed, lastClosed);
        [latb, lonb] = geolinebuf(latClosed, lonClosed, bufwidth, n);
        [lata, lona] = maptrimp(latClosed, lonClosed, latlim, lonlim);
        [lonb, latb] = polybool('intersection', lona, lata, lonb, latb);
    else
        latb = [];
        lonb = [];
    end
else
    if ~isempty(firstClosed)
        % Buffer the edges of closed polygons using a planar topology;
        % convert the input polygons to the same topology and longitude
        % limits; subtract the (interiors of) the closed polygons.
        [latClosed, lonClosed] = subsetCurves( ...
            lat, lon, firstClosed, lastClosed);
        [latc, lonc] = geolinebuf(latClosed, lonClosed, bufwidth, n);
        [lata, lona] = maptrimp(latClosed, lonClosed, latlim, lonlim);
        [lonc, latc] = polybool('subtraction', lonc, latc, lona, lata);
    end
    
    firstOpen = first(~isClosed);
    lastOpen  = last(~isClosed);
    if ~isempty(firstOpen)
        % Buffer the open curves, if any.
        [latOpen, lonOpen] = subsetCurves(lat, lon, firstOpen, lastOpen);
        [lato, lono] = geolinebuf(latOpen, lonOpen, bufwidth, n);
    end
    
    % Combine the results for the open and closed curves. Avoid
    % passing empty inputs to polybool.
    if isempty(firstClosed)
        if isempty(firstOpen)
            latb = [];
            lonb = [];
        else
            latb = lato;
            lonb = lono;
        end
    else
        if isempty(firstOpen)
            latb = latc;
            lonb = lonc;
        else
           [lonb, latb] = polybool('union', lonc, latc, lono, lato);
        end
    end
end

if ~isempty(latb)
    % Convert from planar topology to spherical topology.
    tolSnap = 100*eps(180);
    if any(circumpolar(lonb) ~= 0)
        [lonb, latb] = gluePolygonsOnVerticalEdges( ...
            lonb(:), latb(:), lonlim, tolSnap);
    end
    [latb, lonb] = removeExtraPolarVertices(latb, lonb, tolSnap);
end

% Make shape consistent with input.
latb = latb(:);
lonb = lonb(:);

if rowVectorInput
    latb = latb';
    lonb = lonb';
end

%---------------------------------------------------------------------------

function [x, y] = subsetCurves(x, y, first, last)
% Keep the subset of the NaN-separated curve X,Y defined by the indices
% in first and last. Discard other elements, except for required
% NaN-separators.
    
keep = false(size(x));
for k = 1:numel(first)
    keep(first(k):last(k)) = true;
end
x(~keep) = NaN;
y(~keep) = NaN;
[x, y] = removeExtraNanSeparators(x, y);
