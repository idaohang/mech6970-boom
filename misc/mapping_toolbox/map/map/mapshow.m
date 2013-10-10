function varargout = mapshow(varargin)
%MAPSHOW Display map data without projection
%
%   MAPSHOW(X,Y) or MAPSHOW(X,Y, ..., 'DisplayType', DISPLAYTYPE, ...)
%   displays the map coordinate vectors, X and Y.  X and Y may contain
%   embedded NaNs, delimiting individual lines or polygon parts.
%   DISPLAYTYPE can be 'point', 'line', or 'polygon' and defaults to
%   'line'.
%
%   MAPSHOW(X,Y,Z, ..., 'DisplayType', DISPLAYTYPE, ...) displays a
%   geolocated data grid, where X and Y are M-by-N map coordinate arrays,
%   and Z is an M-by-N array of class double. X, Y, and Z may contain NaN
%   values. DISPLAYTYPE must be set to 'surface', 'mesh', 'texturemap', or
%   'contour'.
%
%   MAPSHOW(Z,R, ..., 'DisplayType', DISPLAYTYPE, ...) displays a regular
%   data grid.  Z is 2-D array of class double and R is a referencing
%   matrix or a spatialref.MapRasterReference object that relates the
%   subscripts of Z to map coordinates. DISPLAYTYPE must be set to
%   'surface', 'mesh', 'texturemap', or 'contour'.  If DISPLAYTYPE is
%   'texturemap', MAPSHOW constructs a surface with ZData values set to 0.
%
%   MAPSHOW(X,Y,I)  
%   MAPSHOW(X,Y,BW) 
%   MAPSHOW(X,Y,A,CMAP) 
%   MAPSHOW(X,Y,RGB)
%   MAPSHOW(..., 'DisplayType', DISPLAYTYPE, ...) displays a geolocated
%   image as a texturemap on a zero-elevation surface.  X and Y are 
%   geolocation arrays in map coordinates and I is a grayscale image, BW is
%   a logical image, A is an indexed image with colormap CMAP, or RGB is a
%   truecolor image. X, Y, and the image array must match in size.  If
%   specified, DISPLAYTYPE must be set to 'image'.  
%
%   MAPSHOW(I,R) 
%   MAPSHOW(BW,R) 
%   MAPSHOW(RGB,R) 
%   MAPSHOW(A,CMAP,R) 
%   MAPSHOW(..., 'DisplayType', DISPLAYTYPE, ...) displays an image
%   georeferenced to map coordinates.  It constructs an image object if the
%   display geometry permits; otherwise, the image is shown as a texturemap
%   on a zero-elevation surface. If specified, DISPLAYTYPE must be set to
%   'image'.
%
%   MAPSHOW(S) or MAPSHOW(S, ..., 'SymbolSpec', SYMSPEC, ...) displays
%   the vector geographic features stored in the geographic data
%   structure S as points, multipoints, lines, or polygons according to
%   the 'Geometry' field of S.  If S is a mapstruct (with 'X' and 'Y'
%   coordinate fields), the coordinate values are used directly to plot
%   features in map coordinates.  If S is a geostruct (with 'Lat' and
%   'Lon' fields), vertices are projected using the Plate Carree
%   projection and a warning is issued.  SYMSPEC specifies the
%   symbolization rules used for the vector data through a structure
%   returned by MAKESYMBOLSPEC.
%
%   MAPSHOW(FILENAME) displays data from FILENAME, according to the type of
%   file format. The DisplayType parameter is set automatically, according
%   to the following table:
%
%       Format                          DisplayType
%       ------                          -----------
%       shapefile                       'point', 'multipoint, 'line', 
%                                       or 'polygon'
%       GeoTIFF                         'image'
%       TIFF/JPEG/PNG with a world file 'image'
%       ARC ASCII GRID                  'surface' (may be overridden)
%       SDTS raster                     'surface' (may be overridden)
%
%   MAPSHOW(AX, ...) sets the axes parent to AX. This is equivalent to 
%   MAPSHOW(..., 'Parent', AX, ...)
%
%   H = MAPSHOW(...) returns a handle to a MATLAB graphics object or, in 
%   the case of polygons, a modified patch object.  If a mapstruct or
%   shapefile name is input, MAPSHOW returns the handle to an hggroup
%   object with one child per feature in the mapstruct or shapefile.  In
%   the case of a polygon mapstruct or shapefile, each child is a modified
%   patch object;  otherwise it is a line object.
%
%   MAPSHOW(..., PARAM1, VAL1, PARAM2, VAL2, ...) specifies parameter/value
%   pairs that modify the type of display or set MATLAB graphics
%   properties. Parameter names can be abbreviated and are case-insensitive.  
%   Refer to the MATLAB Graphics documentation on line, patch, image,
%   surface, mesh, and contour for a complete description of these
%   properties and their values.
%
%   Example 1 
%   ---------
%   % Overlay Boston roads on an orthophoto. 
%   % Includes material (c) GeoEye, all rights reserved.
%   figure
%   mapshow boston.tif
%   axis image off
%
%   % Convert Boston roads to units of survey feet and overlay on
%   % orthophoto.
%   S = shaperead('boston_roads.shp');
%   surveyFeetPerMeter = unitsratio('sf','meter');
%   x = surveyFeetPerMeter * [S.X];
%   y = surveyFeetPerMeter * [S.Y]; 
%   mapshow(x,y)
%
%   Example 2
%   ---------
%   % Display Boston roads and change the LineStyle.
%   roads = shaperead('boston_roads.shp');
%   figure
%   mapshow(roads,'LineStyle',':');
%
%   Example 3 
%   ---------
%   % Display Boston roads using a SymbolSpec.
%   roadspec = makesymbolspec('Line',...
%                             {'ADMIN_TYPE',0,'Color','c'}, ...
%                             {'ADMIN_TYPE',3,'Color','r'},...
%                             {'CLASS',6,'Visible','off'},...
%                             {'CLASS',[1 4],'LineWidth',2});
%   figure
%   mapshow('boston_roads.shp','SymbolSpec',roadspec);
%  
%   Example 4 
%   ---------
%   % Override a graphics property of the SymbolSpec.
%   roadspec = makesymbolspec('Line',...
%                             {'Default', 'Color', 'yellow'}, ...
%                             {'ADMIN_TYPE',0,'Color','c'}, ...
%                             {'ADMIN_TYPE',3,'Color','r'},...
%                             {'CLASS',6,'Visible','off'},...
%                             {'CLASS',[1 4],'LineWidth',2});
%   figure
%   mapshow('boston_roads.shp', 'Color', 'black', 'SymbolSpec', roadspec);
%
%   Example 5 
%   ---------
%   % Override default properties of the line.
%   roadspec = makesymbolspec('Line',...
%                             {'Default', 'Color', 'black'}, ...
%                             {'ADMIN_TYPE',0,'Color','c'}, ...
%                             {'ADMIN_TYPE',3,'Color','r'},...
%                             {'CLASS',6,'Visible','off'},...
%                             {'CLASS',[1 4],'LineWidth',2});
%   figure
%   mapshow('boston_roads.shp','SymbolSpec',roadspec);
%
%   Example 6 
%   ---------
%   % Display an orthophoto including a pond with three large islands.
%   [ortho, cmap] = imread('concord_ortho_w.tif');
%   R = worldfileread('concord_ortho_w.tfw', 'planar', size(ortho));
%   figure
%   mapshow(ortho, cmap, R)
%
%   % Overlay a polygon representing the same pond 
%   % (feature 14 in the concord_hydro_area shapefile).  
%   % Note that the islands are visible in the orthophoto 
%   % through three "holes" in the pond polygon. 
%   pond = shaperead('concord_hydro_area.shp', 'RecordNumbers', 14);
%   mapshow(pond, 'FaceColor', [0.3 0.5 1], 'EdgeColor', 'black')
%
%   % Overlay roads in the same figure.
%   mapshow('concord_roads.shp', 'Color', 'red', 'LineWidth', 1);
%
%   Example 7 
%   ---------
%   % Read and view the Mount Washington SDTS DEM terrain data.
%   [Z, R] = sdtsdemread('9129CATD.DDF');
%
%   % View the Mount Washington terrain data as a mesh.
%   figure
%   mapshow(Z, R, 'DisplayType', 'mesh');
%   demcmap(Z)
%
%   % View the Mount Washington terrain data as a surface.
%   figure
%   mapshow(Z, R, 'DisplayType', 'surface');
%   demcmap(Z)
%
%   % View as a 3-D surface.
%   view(3);
%   axis normal
%
%   Example 8 
%   ----------
%   % Display the grid and contour lines of 
%   % Mount Washington and Mount Dartmouth.
%
%   % Read the terrain data files.
%   [Z_W, R_W] = arcgridread('MtWashington-ft.grd');
%   [Z_D, R_D] = arcgridread('MountDartmouth-ft.grd');
%
%   % Display the terrain data as a surface.
%   figure('Renderer', 'zbuffer')
%   hold on
%   mapshow(Z_W, R_W, 'DisplayType', 'surface');
%   mapshow(Z_D, R_D, 'DisplayType', 'surface');
%
%   % Overlay black contour lines with labels onto the surface.
%   % Set the Z values of the contours to the maximum value of the
%   % corresponding surface.
%   cW = mapshow(Z_W, R_W, 'DisplayType', 'contour', ...
%      'LineColor','black', 'ShowText', 'on');
%   cD = mapshow(Z_D, R_D, 'DisplayType', 'contour', ...
%      'LineColor','black', 'ShowText', 'on');
%   zdatam(get(cW,'Children'), max(Z_W(:)));
%   zdatam(get(cD,'Children'), max(Z_D(:)));
%
%   % Set the colormap appropriate to terrain elevation.
%   demcmap(Z_W)
%
%   See also GEOSHOW, MAKESYMBOLSPEC, MAPVIEW.

% Copyright 1996-2010 The MathWorks, Inc.
% $Revision: 1.1.6.10 $  $Date: 2010/11/17 11:24:25 $

% Verify the number of inputs.
error(nargchk(1,Inf,nargin,'struct'))

% Designate first argument (axes handle) as Parent if present.
varargin = designateAxesArgAsParentArg(mfilename, varargin{:});

% Read data from file if filename provided in argument list and insert
% results at the beginning of argument list.
varargin = importFromFileAndSetDataArgs(mfilename, varargin{:});

% Decide which function to use to display the data (based mainly on the
% data arguments at the beginning of varargin).;
options.showStructFcn = @mapstructshow;
options.showVecFcn    = @mapvecshow;
options.showRasterFcn = @maprastershow;
showFcn = determineShowFcn(options, varargin{:});

% Show the data.
h = showFcn(varargin{:});

% Set the axes appearance properties.
setAxesProperties(h);

% Suppress output if called with no return value and no semicolon.
if nargout == 1
   varargout{1} = h;
end
