function h = globevec(mstruct, lat, lon, height, objectType, varargin) %#ok<INUSL>
%GLOBEVEC Display point, multipoint, line, or polygon on globe
%
%   The HEIGHT input is currently ignored.

% Copyright 2009 The MathWorks, Inc.
% $Revision: 1.1.6.1 $  $Date: 2009/11/09 16:26:05 $

% Verify NaN locations are equal.
assert(isequal(isnan(lat), isnan(lon)), ...
    'map:globevec:inconsistentLatLon', ...
    '%s and %s mismatch in size or NaN locations.','LAT','LON')

if ~isempty(lat) || ~isempty(lon)
    switch(objectType)
        case {'geopoint','geomultipoint'}
            h = globepoint(mstruct, lat, lon, varargin{:});
            
        case 'geoline'
            h = globeline(mstruct, lat, lon, varargin{:});
            
        case 'geopolygon'
            h = globepolygon(mstruct, lat, lon, varargin{:});
    end
else
    % Either xdata or ydata are empty.
    h = reshape([],[0 1]);
end

end

%-----------------------------------------------------------------------

function h = globepoint(mstruct, lat, lon, varargin)
% Construct a line but show only the markers.

[x, y, z] = globe(mstruct, lat(:), lon(:), [], 'geopoint', 'forward');
h = line(x, y, z, 'Marker', '+', ...
    'MarkerEdgeColor', 'red', varargin{:}, 'LineStyle', 'none');
end

%-----------------------------------------------------------------------

function h = globeline(mstruct, lat, lon, varargin)
% Construct a line with a default color of 'blue'.

[x, y, z] = globe(mstruct, lat(:), lon(:), [], 'geoline', 'forward');
h = line(x, y, z, 'Color', 'blue', varargin{:});

end

%-----------------------------------------------------------------------

function h = globepolygon(mstruct, lat, lon, varargin)
% Construct a polygon object for display in a 3-D "globe" system.
% LAT and LON contain the polygon vertices, and use NaNs to
% designate multiple parts, including inner rings.  LAT and LON
% must be the same size and have NaNs in matching locations.  Polygons
% support the same graphics properties as patches.  Return empty if
% LAT and LON are empty.

% Apdated from private/mappolygon.

% Separate out any 'Parent' properties from varargin.
qParent = strncmpi(varargin,'pa',2);
qParent = qParent | circshift(qParent,[0 1]);
parents = varargin(qParent);
varargin(qParent) = [];

if any(~isnan(lat(:))) || any(~isnan(lon(:)))
    [f, vLat, vLon] = geoPolygonToFaceVertex(lat, lon, mstruct.angleunits);
    [vX, vY, vZ] = globe(mstruct, vLat, vLon, [], 'none', 'forward');
    geodata = {'Faces', f, 'Vertices', [vX vY vZ]};
else
    % xdata and ydata both contain nothing but NaN.
    geodata = {'XData', NaN, 'YData', NaN, 'ZData', NaN};
end
% Construct a pale-yellow patch with edges turned off;
% keep it invisible for now.
defaultFaceColor = [1 1 0.5];   % Pale yellow
h = patch(geodata{:}, parents{:}, ...
    'FaceColor',defaultFaceColor, ...
    'EdgeColor','none', ...
    'Visible','off');

% Construct an "edge line" object.
[x, y, z] = globe(mstruct, lat(:), lon(:), [], 'geoline', 'forward');
hEdgeLine = constructEdgeLine(h,x,y,z);

% After setting up listeners and initializing the "update" state, we're
% ready to set the rest of the input properties and make both patch and
% line visible. Also add a delete function to clean up the edge line.
set(h,'Visible','on',varargin{:},'DeleteFcn',@deleteEdgeLine)
set(hEdgeLine,'Visible','on')

%------------------- nested callback function ------------------

    function deleteEdgeLine(hSrc,evnt) %#ok<INUSD>
        % Delete the edge line.
        if ishghandle(hEdgeLine)
            delete(hEdgeLine);
        end
    end

end

%-----------------------------------------------------------------------

function [f, vLat, vLon] = geoPolygonToFaceVertex(lat, lon, angleunits)
% Convert a latitude-longitude polygon to face vertex form. F is the
% face array, vLat is a column vector containing the vertex latitudes,
% and vLon is a column vector, the same size as vLat, containing the
% vertex longitudes.

% Project into a Plate Carree system covering the entire globe,
% triangulate the polygons in the Plate Carree system, and convert the
% vertex locations back to latitude-longitude.
mstruct = defaultm('pcarree');
mstruct = defaultm(mstruct);
[lat, lon] = toDegrees(angleunits, lat, lon);
[x, y] = feval(mstruct.mapprojection, ...
    mstruct, lat(:), lon(:), 'geopolygon', 'forward');
[f, v] = polygonToFaceVertex(x, y);
[vLat, vLon] = feval(mstruct.mapprojection, ...
    mstruct, v(:,1), v(:,2), 'none', 'inverse');
[vLat, vLon] = fromDegrees(angleunits, vLat, vLon);
end
