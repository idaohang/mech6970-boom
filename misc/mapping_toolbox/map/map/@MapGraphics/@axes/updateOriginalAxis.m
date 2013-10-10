function updateOriginalAxis(this)
%UPDATEORIGINALAXIS
%  changes the OrigPosition, OrigXLim and OrigYLim fields of the map axis.

% Copyright 1996-2008 The MathWorks, Inc.
% $Revision: 1.1.6.2 $  $Date: 2008/11/24 14:59:37 $

ax = this.getAxes();

oldunits = get(ax,'Units');
set(ax,'Units','pixels')
this.OrigPosition = get(ax,'Position');
this.OrigXLim = get(ax,'XLim');
this.OrigYLim = get(ax,'YLim');
set(ax,'Units',oldunits)
