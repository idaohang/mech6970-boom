function resizeLimits(this)
% Resize the axes limits to negate any rescaling

% Copyright 2008 The MathWorks, Inc.
% $Revision: 1.1.6.2 $  $Date: 2008/11/24 14:59:32 $

ax = this.getAxes();

p = getpixelposition(ax);
xChngFactor = p(3) / this.OrigPosition(3);
yChngFactor = p(4) / this.OrigPosition(4);

xlim = sum(this.OrigXLim)/2 + [-1/2, 1/2] * diff(this.OrigXLim) * xChngFactor;
ylim = sum(this.OrigYLim)/2 + [-1/2, 1/2] * diff(this.OrigYLim) * yChngFactor;

if diff(xlim) > 0
  set(ax,'XLim',xlim);
end
if diff(ylim) > 0
  set(ax,'YLim',ylim);
end
