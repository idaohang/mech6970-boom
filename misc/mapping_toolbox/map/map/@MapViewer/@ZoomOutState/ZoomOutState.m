function this = ZoomOutState(viewer,mode)

% Copyright 1996-2008 The MathWorks, Inc.
% $Revision: 1.1.6.8 $  $Date: 2008/10/26 14:26:49 $

this = MapViewer.ZoomOutState;

this.MapViewer = viewer;

viewInfo = this.MapViewer.Axis.ViewInfo;
if (isempty(viewInfo))
   this.MapViewer.Axis.ViewInfo = viewInfo;
end

% Turn on zoom but preserve the WindowButtonMotionFcn
WindowButtonMotionFcn = get(viewer.Figure,'WindowButtonMotionFcn');
zoom(viewer.Figure,'outmode');
set(viewer.Figure,'WindowButtonMotionFcn',WindowButtonMotionFcn);            

glassMinusPtr = setptr('glassminus');
viewer.setCursor(glassMinusPtr);
