function renderLayer(this,ax,layername)
%RENDERLAYER Render one layer in the model
%
%   H = RENDERLAYER(AX,LAYERNAME) renders the layer LAYERNAME into the axes
%   AX. 

% Copyright 1996-2007 The MathWorks, Inc.
% $Revision: 1.1.6.3 $  $Date: 2007/11/09 20:26:05 $

for i=1:length(this.Layers)
  names{i} = this.Layers(i).getLayerName;
end
I = strmatch(layername,names,'exact');
if isempty(I)
  error(['map:' mfilename ':mapError'], ...
      'A layer named %s does not exist in this model.',layername)
end
this.Layers(I).render(ax);
