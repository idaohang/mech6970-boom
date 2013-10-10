function scale = getScale(this)
%GETSCALE Returns the absolute map scale.  

% Copyright 1996-2008 The MathWorks, Inc.
% $Revision: 1.1.6.3 $  $Date: 2008/11/24 14:59:27 $

if this.MapUnitInCM == 0
    scale = [];
else
    ax = this.getAxes();
    units = get(ax,'Units');
    set(ax,'Units','Centimeters')
    p = get(ax,'Position');
    set(ax,'Units',units)
    xScale = p(3) / (diff(get(ax,'XLim')) * this.MapUnitInCM);
    yScale = p(4) / (diff(get(ax,'YLim')) * this.MapUnitInCM);
    scale = min([xScale, yScale]);
end
