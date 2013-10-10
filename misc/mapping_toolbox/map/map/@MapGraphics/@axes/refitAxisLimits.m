function refitAxisLimits(this)
%REFITAXISLIMITS
%   sets the axis limits so that the a 10 point border is maintained
%   around the axis. 

% Copyright 1996-2008 The MathWorks, Inc.
% $Revision: 1.1.6.2 $  $Date: 2008/11/24 14:59:31 $

ax = this.getAxes();

p = getpixelposition(ax);
if prod(p(3:4)) > 0
    axisProportion = p(3)/p(4);
end

oldXlim = get(ax,'XLim');
oldYlim = get(ax,'YLim');
if (diff(oldXlim)/diff(oldYlim) < 0)
    return;
end

if ( diff(oldXlim)/diff(oldYlim) < axisProportion )
  % the x Axis lims need rescaling
  newLimRange = diff(oldYlim) * p(3)/p(4);
  set(ax,'XLim', sum(oldXlim)/2 + newLimRange * [-1/2, 1/2]);
else
  newLimRange = diff(oldXlim) * p(4)/p(3);
  set(ax,'YLim', sum(oldYlim)/2 + newLimRange * [-1/2, 1/2]);
end
