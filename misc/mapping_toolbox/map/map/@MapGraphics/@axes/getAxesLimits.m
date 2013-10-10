function bbox = getAxesLimits(this)
%GETAXESLIMITS
%
%   Returns the axes lower-left and upper-right corners [lower-left-x,y;
%   upper-right-x,y],
%   
%      or equivalently,  [left      bottom;
%                         right        top]

% Copyright 1996-2008 The MathWorks, Inc.
% $Revision: 1.1.6.3 $  $Date: 2008/11/24 14:59:24 $

ax = this.getAxes();
xLimits = get(ax,'XLim');
yLimits = get(ax,'YLim');
bbox = [xLimits(:) yLimits(:)];
