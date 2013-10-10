function zoom(this,varargin)
%ZOOM   Zoom in and out on a map.
%   This method simply acts as a wrapper to call the HG zoom tool
%   and it ensures that all the layers are considered for when 
%   determining the default/reset plot view.

% Copyright 1996-2008 The MathWorks, Inc. 
% $Revision: 1.1.6.4 $  $Date: 2008/11/24 15:00:11 $

viewInfo = this.Axis.ViewInfo;
if (isempty(viewInfo))
    this.Axis.ViewInfo = localCreateViewInfo(this);
end

% Call the HG zoom
zoom(varargin{:});

function [viewinfo] = localCreateViewInfo(this)  
% Same as the resetplotview localCreateViewInfo except that
% the XLim and YLim are from the largest layers on the map
% axis

% Turn off the listeners 
scribefiglisten(this.Figure,'off');

% Create a temporary Map axis for obtaining the maximum
% axis limits for the layers currently on the map.
hAxes = this.getAxes();
tmpAxes = MapGraphics.axes('Parent',get(hAxes,'Parent'),'Visible','off');
tmpAxes.setAxesLimits(this.map.getBoundingBox.getBoxCorners);
tmpAxes.refitAxisLimits(); 

viewinfo.XLim = get(tmpAxes.getAxes(),'XLim');
viewinfo.YLim = get(tmpAxes.getAxes(),'YLim');

% Store axes view state
viewinfo.DataAspectRatio = get(hAxes,'DataAspectRatio');
viewinfo.DataAspectRatioMode = get(hAxes,'DataAspectRatioMode');
viewinfo.PlotBoxAspectRatio = get(hAxes,'PlotBoxAspectRatio');
viewinfo.PlotBoxAspectRatioMode = get(hAxes,'PlotBoxAspectRatioMode');
viewinfo.XLimMode = get(hAxes,'XLimMode');
viewinfo.YLimMode = get(hAxes,'YLimMode');
viewinfo.ZLim = get(hAxes,'zLim');
viewinfo.ZLimMode = get(hAxes,'ZLimMode');
viewinfo.CameraPosition = get(hAxes,'CameraPosition');
viewinfo.CameraViewAngleMode = get(hAxes,'CameraViewAngleMode');
viewinfo.CameraTarget = get(hAxes,'CameraTarget');
viewinfo.CameraPositionMode = get(hAxes,'CameraPositionMode');
viewinfo.CameraUpVector = get(hAxes,'CameraUpVector');
viewinfo.CameraTargetMode = get(hAxes,'CameraTargetMode');
viewinfo.CameraViewAngle = get(hAxes,'CameraViewAngle');
viewinfo.CameraUpVectorMode = get(hAxes,'CameraUpVectorMode');

delete(tmpAxes.getAxes())
% Turn listeners back on 
scribefiglisten(this.Figure,'on');
