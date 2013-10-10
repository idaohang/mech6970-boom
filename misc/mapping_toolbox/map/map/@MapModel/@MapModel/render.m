function render(this,ax)
%RENDER Render all layers in the model
%
%   RENDER(AX) renders all the layers in the model into the axes AX
%   according to the order defined in the configuration of the model.
%

%   Copyright 1996-2006 The MathWorks, Inc.
%   $Revision: 1.1.6.3 $  $Date: 2006/12/10 20:04:28 $

if ~isempty(this.Configuration)
  unorderedLayers = get(this,'layers');

  % Don't render empty layers
  for i=1:length(unorderedLayers)
    if isempty(unorderedLayers(i))
      unorderedLayers(i) = [];
    end
  end

  % Render layers in order of the Configuration
  if length(unorderedLayers) == 1
    idx = 1;
  else
    [sortedUnorderedNames,idx] = sort(get(unorderedLayers,'layername'));
  end
  sortedUnordered = unorderedLayers(idx);
  [sortedOrderedNames,idx] = sort(this.Configuration);
  % The FLIPUD is because the layer orders are stored topmost first, but the
  % topmost layer must be rendered last.
  sortedLayers = sortedUnordered(flipud(idx(:)));

  % Let the layers render themselves
  for i=1:length(sortedLayers)
    sortedLayers(i).render(ax);
  end
end