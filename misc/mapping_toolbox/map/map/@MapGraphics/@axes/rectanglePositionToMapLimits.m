function [xLimits, yLimits] = rectanglePositionToMapLimits(this,rect_pos)
% Given a rectangle defined by a position vector in pixel units, compute
% the corresponding limits in map coordinates.

% Copyright 2008 The MathWorks, Inc.
% $Revision: 1.1.6.2 $  $Date: 2008/11/24 14:59:30 $

ax = this.getAxes();

xLimits = get(ax,'XLim');
yLimits = get(ax,'YLim');

axes_pos = getpixelposition(ax);

xlim_pos_ratio = diff(xLimits)/(axes_pos(3));
ylim_pos_ratio = diff(yLimits)/(axes_pos(4));

xLimits = xLimits(1) + xlim_pos_ratio * (rect_pos(1) + [0 rect_pos(3)]);
yLimits = yLimits(1) + ylim_pos_ratio * (rect_pos(2) + [0 rect_pos(4)]);
