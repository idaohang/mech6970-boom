function setLayerLegend(this,layername,legend,axis)
%SETLAYERLEGEND Sets a legend for a specific layer
%
%   SETLAYERLEGEND(LAYERNAME,LEGEND,AXIS) Sets LAYERNAME's legend to the
%   mapmodel.legend LEGEND, updating the corresponding properties of the axis
%   (MapGraphics.axes) children.

%   Copyright 1996-2009 The MathWorks, Inc.
%   $Revision: 1.1.6.4 $  $Date: 2009/11/09 16:25:50 $

layer = this.getLayer(layername);
layer.setLegend(legend);

handles = axis.getLayerHandles(layer.getLayerName);
for i=1:length(handles)
  % Make a "fake" feature using the real feature's attributes
  tmpfeature.Attributes = getAttributes(handles(i));
  props = legend.getGraphicsProperties(tmpfeature);
  set(handles(i),props)
end

function attributes = getAttributes(h)
attributes = getappdata(h,'Attributes');

