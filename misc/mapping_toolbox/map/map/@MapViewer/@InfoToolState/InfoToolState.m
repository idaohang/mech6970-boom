function this = InfoToolState(viewer)

%   Copyright 1996-2007 The MathWorks, Inc.
%   $Revision: 1.1.6.3 $  $Date: 2007/10/10 20:47:56 $


if isempty(viewer.PreviousInfoToolState)
  this = MapViewer.InfoToolState;
else
  this = viewer.PreviousInfoToolState;
  viewer.PreviousInfoToolState = [];
end

this.Viewer = viewer;

viewer.setCursor({'Pointer','crosshair'});
                  
this.setActiveLayer(viewer, viewer.ActiveLayerName);