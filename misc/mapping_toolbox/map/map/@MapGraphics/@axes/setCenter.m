function setCenter(this,center)
%SETCENTER Set axis center
%
%   SETCENTER(CENTER) sets the center of the map to be CENTER. CENTER is a 2
%   element array [X Y], the x and y coordinates, in map units, of the center of
%   the axes.

% Copyright 1996-2008 The MathWorks, Inc.
% $Revision: 1.1.6.4 $  $Date: 2008/11/24 14:59:35 $

if isempty(this.mapUnitInCM)
    error(['map:' mfilename ':mapError'], ...
        'MapUnitInCM must be set before changing the center.')
else
    ax = this.getAxes();
    units = get(ax,'Units');
    set(ax,'Units','Centimeters')
    p = get(ax,'Position');
    xLimits = center(1) + p(3) * [-1/2 1/2] / (this.getScale * this.mapUnitInCM);
    yLimits = center(2) + p(4) * [-1/2 1/2] / (this.getScale * this.mapUnitInCM);
    set(ax,'XLim',xLimits,'YLim',yLimits)
    set(ax,'Units',units)
end
