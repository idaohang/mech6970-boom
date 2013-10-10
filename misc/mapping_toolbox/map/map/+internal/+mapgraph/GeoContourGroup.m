%GeoContourGroup Contours in latitude-longitude
%
%       FOR INTERNAL USE ONLY -- This class is intentionally
%       undocumented and is intended for use only within other toolbox
%       classes and functions. Its behavior may change, or the class
%       itself may be removed in a future release.
%
%   GeoContourGroup properties:
%      ContourLabels - Whether to label contours and which to label
%      Fill - Color areas between contour lines
%      FillAlpha - Transparency of contour-fill polygons
%      FillColor - Value or method for selecting contour-fill polygon colors
%      FillColormap - Color map for filled contour intervals
%      FillZ - Height at which to display contour-fill polygons
%      LabelSpacing - Distance between labels in points
%      LevelList - Vector of levels at which contours are computed
%      LineColor - Color value or method for selecting contour line colors
%      LineColormap - Color map for contour lines
%      LineStyle - LineStyle property for contour lines
%      LineWidth - Width of contour lines in points
%      LineZ - Height at which to display contour lines
%      SpatialRef - Geographic spatial referencing object or structure
%      ZData - Data grid from which contour lines are generated
%
%   GeoContourGroup methods:
%      GeoContourGroup - Construct GeoContourGroup object
%      delete - Delete GeoContourGroup object
%      fillPolygonColors - Return one color per contour interval
%      contourLineColors - Return one color per contour level
%      getContourLines - Geostruct array with one line per level
%      getFillPolygons - Geostruct array with one polygon per contour interval
%      getTextLabelHandles - Find associated text objects
%      refresh - Create or update contour display
%      reproject - Refresh display in response to map axes geometry

% Copyright 2010 The MathWorks, Inc.
% $Revision: 1.1.6.5 $  $Date: 2010/09/24 14:33:38 $

classdef GeoContourGroup < internal.mapgraph.ContourGroup
        
    %------------------------ Public methods ---------------------------
    
    methods
        
        function h = GeoContourGroup(ax, Z, R, levels)
            %GeoContourGroup Construct GeoContourGroup object 
            %
            %   h = GeoContourGroup(ax, Z, R, levels) constructs a
            %   GeoContourGroup object h, and an associated hggroup object,
            %   given a parent axes AX, data grid Z, referencing object R,
            %   and a list of contour levels.
            
            % Use base class to initialize new object.
            h = h@internal.mapgraph.ContourGroup(ax, Z, R, levels);                        
        end

        
        function L = getContourLines(h)
            %getContourLines Contour lines in latitude-longitude
            %
            %   L = getContourLines(h) returns a line geostruct L with
            %   contour lines corresponding to the current ZData,
            %   SpatialRef, and LevelList properties of the contour
            %   object h, in a geographic coordinate system. There is an
            %   element for each contour level intersected by the range
            %   of values in h.ZData, and the contour level values are
            %   stored in a 'Level' field.
            
            if isempty(h.pContourLines)
                % For efficiency, call geocontours only if the contour
                % lines structure hasn't already been computed.
                h.pContourLines ...
                    = geocontours(h.ZData, h.SpatialRef, h.LevelList);
            end
            L = h.pContourLines;
        end
        
        
        function P = getFillPolygons(h)
            %getFillPolygons Fill polygons in latitude-longitude
            %
            %   P = getFillPolygons(h) returns a polygon geostruct P
            %   with fill polygons corresponding to the current ZData,
            %   SpatialRef, and LevelList properties of the contour
            %   object h, in a geographic coordinate system. There is
            %   one element for each contour interval intersected by the
            %   range of values in h.ZData, and the limits of each
            %   contour interval are stored in 'MinLevel' and 'MaxLevel'
            %   fields.

            if isempty(h.pFillPolygons)
                % For efficiency, call geocontours only if the fill
                % polygons structure hasn't already been computed.
               [h.pContourLines, h.pFillPolygons] ...
                    = geocontours(h.ZData, h.SpatialRef, h.LevelList);
            end
            P = h.pFillPolygons;
        end
        
        function reproject(h)
            %reproject Reproject contour lines and fill
            %
            %   reproject(h) Refreshes the display in response to changes
            %   in the geometric properties of the map axes ancestor of the
            %   hggroup associated with the contour object h.
            
            h.refresh()
        end
        
    end
    
    %------------------ Private and protected methods ---------------------
  
    methods (Access = protected)
                
        function hLine = constructContourLine(h, S, zdata, varargin)
            hLine = projectLine(h.HGGroup, S.Lat, S.Lon, zdata, varargin{:});
        end
        
        
        function hPolygon = constructFillPolygon(h, S, zdata, varargin)
            hPolygon = projectPolygonFaces(h.HGGroup, S.Lat, S.Lon, zdata, varargin{:});
        end
        
        
        function validateOnRefresh(h)
            % Check this also: isequal(size(Z), h.SpatialRef.RasterSize)
            if h.usingGlobe()
                if ~strcmpi(h.ContourLabels,'none')
                    warning('map:GeoContourGroup:labelsWithGlobe', ...
                        'Contour labels will not display in a ''%s'' axes.', ...
                        'globe')
                end
                if strcmpi(h.LineZ, 'levels')
                    warning('map:GeoContourGroup:elevateWithGlobe', ...
                        ['Contours might not display properly in 3-D' ...
                        ' in a ''%s'' axes.'], 'globe')
                end
            end
        end
        
        
        function tf = labelsPermitted(h)
            tf = ~h.usingGlobe();
        end
        
        
        function validateSpatialRef(~, value)
            validateattributes(value, {'double','struct', ...
                'spatialref.GeoRasterReference'},{'nonempty'})
        end
        
    end
    
    methods (Access = private)
        
        function tf = usingGlobe(h)
            % True if and only if the axes ancestor of the hggroup is a
            % map axes with MapProjection set to 'globe'.
            ax = ancestor(h.HGGroup,'axes');
            tf= ismap(ax) && strcmp(getm(ax,'MapProjection'),'globe');
        end
        
    end
    
end

%-------------------------- Non-Method Functions -------------------------

function h = projectLine(parent, lat, lon, zdata, varargin)
% Project and display line objects.

ax = ancestor(parent,'axes');
if ismap(ax)
    % There is a map axes; trim and project.
    mstruct = gcm(ax);
    usingGlobe = strcmp(mstruct.mapprojection,'globe');
    if ~usingGlobe
        [x, y] = feval(mstruct.mapprojection, mstruct, ...
            lat, lon, 'geoline', 'forward');
        z = zdata + zeros(size(x));
    else
        [x, y, z] = globe(mstruct, ...
            lat, lon, zdata + zeros(size(lat)), 'geoline', 'forward');
    end
else
    % No map axes; trim with maptriml then display in ordinary axes.
    [y, x] = maptriml(lat, lon, [-90 90], [-180 180]);
    
    % Ensure row vectors.
    x = x(:)';
    y = y(:)';
    z = zdata + zeros(size(x));
end

% Construct the (multipart) contour line.
h = line('XData', x, 'YData', y, 'ZData', z, varargin{:}, 'Parent', parent);

end

%--------------------------------------------------------------------------

function h = projectPolygonFaces(parent, lat, lon, zdata, varargin)
% Project and display polygon faces using a patch object

f = [];

ax = ancestor(parent,'axes');
if ismap(ax)
    % There is a map axes; trim and project.
    mstruct = gcm(ax);
    usingGlobe = strcmp(mstruct.mapprojection,'globe');
    if ~usingGlobe
        [x, y] = feval(mstruct.mapprojection, mstruct, ...
            lat, lon, 'geopolygon', 'forward');
        if any(~isnan(x))
            [f,v] = polygonToFaceVertex(x,y);
            if ~isempty(f)
                v(:,3) = zdata;
            end
        end
    else
        [f, vLat, vLon] = geoPolygonToFaceVertex(lat, lon);
        [vX, vY, vZ] = globe(mstruct, vLat, vLon, zdata, 'none', 'forward');
        v = [vX vY vZ];
    end
else
    % No map axes; use maptrimp to convert to planar topology then display.
    [y, x] = maptrimp(lat, lon, [-90 90], [-180 180]);
    if any(~isnan(x))
        [f,v] = polygonToFaceVertex(x,y);
    end
end

if ~isempty(f)
    % Construct a patch without visible edges.
    h = patch('Faces', f, 'Vertices', v, varargin{:}, 'Parent', parent);
else
    h = [];
end

end

%-----------------------------------------------------------------------

function [f, vLat, vLon] = geoPolygonToFaceVertex(lat, lon)
% Convert a latitude-longitude polygon to face vertex form. F is the
% face array, vLat is a column vector containing the vertex latitudes,
% and vLon is a column vector, the same size as vLat, containing the
% vertex longitudes.
%
%     Adapted from subfunction in toolbox/map/map/private/globevec

% Project into a Plate Carree system covering the entire globe,
% triangulate the polygons in the Plate Carree system, and convert the
% vertex locations back to latitude-longitude.
mstruct = defaultm('pcarree');
mstruct = defaultm(mstruct);
[x, y] = feval(mstruct.mapprojection, ...
    mstruct, lat(:), lon(:), 'geopolygon', 'forward');
[f, v] = polygonToFaceVertex(x, y);
if ~isempty(f)
    [vLat, vLon] = feval(mstruct.mapprojection, ...
        mstruct, v(:,1), v(:,2), 'none', 'inverse');
else
    vLat = [];
    vLon = [];
end

end

%--------------------------------------------------------------------------

function [faces, vertices] = polygonToFaceVertex(x,y)
% Face-vertex decomposition of multipart polygon based on constrained
% Delaunay triangulation.
%
%    Replicated from toolbox/map/map/private/polygonToFaceVertex

if ~isempty(x)
    % Construct a constrained Delaunay triangulation.
    dt = internal.map.polygonToDelaunayTri(x, y);
    
    % Get the vertices and faces.
    vertices = dt.X;
    faces = dt.Triangulation;
    
    % Discard faces that fall outside the polygon.
    inside = dt.inOutStatus();
    faces(~inside,:) = [];
else
    faces = [];
    vertices = [];
end

end
