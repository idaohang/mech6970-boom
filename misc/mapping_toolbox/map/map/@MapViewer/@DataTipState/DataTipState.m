function this = DataTipState(viewer)

% Copyright 1996-2008 The MathWorks, Inc.
% $Revision: 1.1.6.7 $  $Date: 2008/10/02 18:56:27 $

if ~isempty(viewer.PreviousDataTipState)
    this = viewer.PreviousDataTipState;
    viewer.PreviousDataTipState = [];
else
    this = MapViewer.DataTipState;
end

this.Viewer = viewer;
                  
viewer.setCursor({'Pointer','crosshair'});

this.ActiveLayerDisplay = viewer.DisplayPane.ActiveLayerDisplay;

this.setActiveLayer(viewer,viewer.ActiveLayerName);
