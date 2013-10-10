function setScale(this,scale)
%SETSCALE set map scale
%
%   SETSCALE(DISPLAYSCALE) sets the absolute map scale to be the value of
%   DISPLAYSCALE. Typically this number is small. For example, setting the
%   DISPLAYSCALE to be 0.00001 means that 1 mile on the ground covers 0.0001
%   miles on the map.

% Copyright 1996-2008 The MathWorks, Inc.
% $Revision: 1.1.6.5 $  $Date: 2008/11/24 14:59:36 $

if isempty(this.MapUnitInCM)
    error(['map:' mfilename ':mapError'], ...
        'MapUnitInCM must be set before changing the scale.')
else
    ax = this.getAxes();
    units = get(ax,'Units');
    set(ax,'Units','Centimeters')    
    xCenterInMapUnits = sum(get(ax,'XLim')) / 2;
    yCenterInMapUnits = sum(get(ax,'YLim')) / 2;
    p = get(ax,'Position');    
    xLimits = xCenterInMapUnits + p(3) * [-1/2 1/2] / (scale * this.MapUnitInCM);
    yLimits = yCenterInMapUnits + p(4) * [-1/2 1/2] / (scale * this.MapUnitInCM);
    set(ax,'XLim',xLimits,'YLim',yLimits)
    set(ax,'Units',units)
    
    this.updateOriginalAxis();
end
