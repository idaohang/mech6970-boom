function render(this,layerName,legend,ax,visibility)
%RENDER Render the polygon component.
%
%   RENDER(LAYERNAME, LEGEND, AX, VISIBILITY) renders all features
%   of the polygon component into the axes AX using the symbolization
%   defined in the legend, LEGEND, for the layer defined by LAYERNAME.
%   The polygon visibility is defined by VISIBILITY.

% Copyright 1996-2008 The MathWorks, Inc.
% $Revision: 1.1.6.8 $  $Date: 2008/11/24 14:59:46 $

mappolygon = mapgate('mappolygon');

features = this.Features;
for k = 1:numel(features)
    % Polygon vertex arrays.
    xdata = features(k).xdata;
    ydata = features(k).ydata;
    
    % Graphics properties from symbolization rules.
    properties = legend.getGraphicsProperties(features(k));
    
    % Construct the k-th polygon --
    % a patch with an associated "edge line".
    h = mappolygon(xdata, ydata, ...
        'Tag', layerName, ...
        'Parent', ax, ...
        'Visible', visibility, ...
        'HitTest', 'off', ...
        properties);
    
    % Set 'HitTest' on the edge line.
    hEdgeLine = getappdata(h,'EdgeLine');
    set(hEdgeLine,'HitTest','off')
    
    % Store the Attributes structure in the appdata of the patch.
    setappdata(h,'Attributes',features(k).Attributes)
end
