function h = getLayerHandles(this,layerName)
%Get all graphics handles for a layer with the specified layerName.

% Copyright 1996-2008 The MathWorks, Inc.
% $Revision: 1.1.6.6 $  $Date: 2008/11/24 14:59:26 $

% Group results by name, and be sure to return handles (not doubles).
children = get(this.getAxes(),'Children');
h = children(strmatch(layerName, ...
        arrayfun(@getLayerName,children,'UniformOutput',false),'exact'));

%-----------------------------------------------------------------------

function layerName = getLayerName(hChild)
% Get the layer name, if any, associated with a given child object.

layerName = get(hChild,'Tag');
