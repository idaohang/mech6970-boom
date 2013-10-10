function h = renderMapData(ax, mapdata)
%RENDERMAPDATA Render the mapdata structure onto the axes. 
%
%   H = RENDERMAPDATA(AX, MAPDATA) renders the MAPDATA structure onto the
%   axes AX.  MAPDATA will hold the rendering function name as well as the
%   function arguments. For a complete description of the MAPDATA structure
%   see PARSEMAPINPUTS.
%  
%   See also BUILDMAPDATA, MAPSHOW, MAPVIEW, PARSEMAPINPUTS, READMAPDATA.

%   Copyright 1996-2005 The MathWorks, Inc.
%   $Revision: 1.1.6.3 $  $Date: 2005/06/20 03:10:54 $

%-------------------------------------------------------------------------

% Render the map's data onto the axes.

if isempty(mapdata.renderFcn)

  switch mapdata.renderFcnName

     case 'surface'
       mapdata.renderFcn = @surface;

     case 'surfacem'
       mapdata.renderFcn = @surfacem;

     case 'mesh'
       mapdata.renderFcn = @mesh;

     case 'meshm'
       mapdata.renderFcn = @meshm;

     case 'renderRegularGeoDataGrid'
       mapdata.renderFcn = @renderRegularGeoDataGrid;

     case 'renderComponent'
       mapdata.renderFcn = @renderComponent;

     case 'renderImageComponent'
       mapdata.renderFcn = @renderImageComponent;

     case 'contourwrap'
       mapdata.renderFcn = @contourwrap;

     otherwise
       % This line should never execute because all
       % possible values of renderFcnName have already been enumerated above.
       mapdata.renderFcn = str2func(mapdata.renderFcnName);
   end
end
h = mapdata.renderFcn(mapdata.args{:});

% Set the Handle Graphics parameter value pairs,
%  including the Parent
if isfield(mapdata,'HGparams') && ~isempty(mapdata.HGparams)
  set(h,mapdata.HGparams);
end

% Set the axes to equal or image
setMapAxes(ax, mapdata);

%----------------------------------------------------------------------
function setMapAxes(ax, mapdata)
% Set the axes depending on the type of data being rendered.

% If hold is on, do nothing
if ishold
  return
end

% Set the map's axes according to its type.
switch lower(mapdata.type)
   
   case  'image'
      xLimMode = get(gca, 'XLimMode');
      yLimMode = get(gca, 'YLimMode');
      axis(ax,'image');
      if ~isequal(xLimMode, get(gca,'XLimMode') ) || ...
         ~isequal(yLimMode, get(gca,'YLimMode') )
         % The mode for the limits has changed
         % Reset back to the original mode
         set(gca,'XLimMode', xLimMode, 'YLimMode', yLimMode);
       end
     

   case 'grid'
      view(2);
      if ~strcmp(mapdata.renderFcnName(end),'m')
       % Do not set for the 'm'  ending functions
         if strcmp(mapdata.renderFcnName,'texturemap')
            axis(ax,'equal','fill') ;
         else
            axis(ax,'equal') ;
         end
      end

   otherwise
      % Set the axes equal for vector data
      axis(ax,'equal') ;
end

%----------------------------------------------------------------------
function h = renderComponent(component, layerName,rules,ax,vis)
% Renders a MapGraphics component onto the axes ax.

h = component.render(layerName,rules,ax,vis);

%----------------------------------------------------------------------
function h = contourwrap(varargin)
% Wraps the contour function to return the handle h.

fcn = varargin{1};
[c,h] = fcn(varargin{2:end});

%----------------------------------------------------------------------
function h = renderImageComponent(I, XData, YData)
% Render an image using the Handle Graphics image function.

h = image('CData',I,'XData',XData,'YData',YData);

%----------------------------------------------------------------------
function h = renderRegularGeoDataGrid(R, x, y, Z);
% Renders a geographic regular data grid onto the axes ax.

h = surfacem(x, y, Z);
userData.maplegend = R;
set(h, 'UserData', userData);

