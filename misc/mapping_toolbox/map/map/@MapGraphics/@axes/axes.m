function this = axes(varargin)
%AXES 
%
%  AXES is a subclass of the HG Axes object. The axes is a viewer for the
%  mapmodel.mapmodel object, MODEL.

% Copyright 1996-2008 The MathWorks, Inc.
% $Revision: 1.1.6.3 $  $Date: 2008/11/24 14:59:22 $

this = MapGraphics.axes;

this.hAxes = axes(varargin{:});

setappdata(this.getAxes(),'MapGraphicsAxesObject',this)

set(this.getAxes(),...
    'DataAspectRatioMode','manual',...
    'DataAspectRatio',[1 1 1],...
    'PlotBoxAspectRatioMode','auto',...
    'NextPlot','Add',...
    'ALimMode','manual',...
    'CLimMode','manual',...
    'XLimMode','manual',...
    'YLimMode','manual',...
    'ZLimMode','manual',...
    'DeleteFcn', @deleteAxes)

    %---------------------- Nested callback function ----------------
    
    function deleteAxes(hSrc,evnt)  %#ok<INUSD>
        if ishandle(this) % Applying ishandle to a non-HG object
            delete(this)
        end
    end
end
