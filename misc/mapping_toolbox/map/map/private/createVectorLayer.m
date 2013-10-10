function [layer,component] = createVectorLayer(S, name)
%CREATEVECTORLAYER Create a vector layer from a geographic data structure.
%
%   [LAYER, COMPONENT] = CREATEVECTORLAYER(S, NAME) creates a layer, LAYER,
%   and a component, COMPONENT, given a geographic data structure, S, and
%   the layer name, NAME.

% Copyright 1996-2008 The MathWorks, Inc.
% $Revision: 1.1.6.3 $  $Date: 2008/11/24 15:00:16 $

attributes = geoattribstruct(S);
if ~isempty(attributes)
  attrnames = fieldnames(attributes);

  % add Index to the geographic data struct as the 
  % the last attribute, if it does not exist.
  if isempty(strmatch('index',lower(attrnames),'exact'))
     attrnames = [attrnames; {'INDEX'}];
     indexNums = num2cell(1:length(attributes));
     [attributes(1:length(indexNums)).INDEX] = deal(indexNums{:});
  end
else
  attrnames = '';
end

% Assume only one Geometry type per structure.
type = S(1).Geometry;
switch lower(type)
  case {'point', 'multipoint'}
     layer = MapModel.PointLayer(name);
     component = MapModel.PointComponent(attrnames);
  case 'line'
     layer = MapModel.LineLayer(name);
     component = MapModel.LineComponent(attrnames);
  case 'polygon'
     layer = MapModel.PolygonLayer(name);
     component = MapModel.PolygonComponent(attrnames);
  otherwise
     error(['map:' mfilename ':invalidGeometryType'], ...
         'Geometry type ''%s'' is not supported.',type)
end
component.addFeatures(S,attributes);
layer.addComponent(component);
