function center = getCenter(this)
%GETCENTER Get center of axes
%
%   CENTER = GETCENTER returns the coordinates of the center of the map in map
%   units. CENTER is a 2 element array [X Y], the x and y coordinates of the
%   center of the axes.

% Copyright 1996-2008 The MathWorks, Inc.
% $Revision: 1.1.6.3 $  $Date: 2008/11/24 14:59:25 $

ax = this.getAxes();
center = [sum(get(ax,'XLim')) sum(get(ax,'YLim'))]/2;
