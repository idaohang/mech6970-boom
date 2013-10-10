function localPan(this,startPt)
% Shift axes limits by the difference between the current point and an
% input startPt.

% Copyright 2008 The MathWorks, Inc.
% $Revision: 1.1.6.2 $  $Date: 2008/11/24 14:59:29 $

ax = this.getAxes();
p = get(ax,'CurrentPoint');
set(ax,...
    'XLim',get(ax,'XLim') - (p(1) - startPt(1)), ...
    'YLim',get(ax,'YLim') - (p(3) - startPt(2)));  
