function adjustPositionInPoints(this, positionDelta)
% Use 4-by-1 positionDelta vector to adjust the axis position in points.

% Copyright 1996-2008 The MathWorks, Inc.
% $Revision: 1.1.6.2 $  $Date: 2008/11/24 14:59:21 $

ax = this.getAxes();
oldunits = get(ax,'Units');
set(ax,'Units','Points')
set(ax,'Position',get(ax,'Position') + positionDelta)
set(ax,'Units',oldunits)
