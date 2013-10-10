function setAxesLimits(this,bbox)
%SETAXESLIMITS Set extent of axes
%
%   SETAXESLIMITS(BBOX) sets the extent of the axes to the bounding box BBOX.

% Copyright 1996-2008 The MathWorks, Inc.
% $Revision: 1.1.6.4 $  $Date: 2008/11/24 14:59:34 $

x = squeeze(bbox(:,1,:)); % 2-by-n with all the x-bounds
y = squeeze(bbox(:,2,:)); % 2-by-n with all the y-bounds

set(this.getAxes(),...
    'XLim',roundinterval([min(x(:)) max(x(:))]),...
    'YLim',roundinterval([min(y(:)) max(y(:))]));

%--------------------------------------------------------------------------------

function y = roundinterval( x )

% Round out an interval: Let d be the length of the interval divided by 10.
% Find f, the closest value to d of the form 10^n, 2 * 10^n, or 5 * 10^n.
% Subtract f from min(x) and round down to the nearest multiple of f.
% Add f to max(x) and round up to the nearest multiple of f.

d = abs(x(2)-x(1))/10;
if d ==0
    % If the interval x degrades to a point, numerically nudge the interval
    % x to 2*eps(x) about the location of the point. HG requires that XLim
    % and YLim intervals are increasing.
    point_val = x(1);
    y = [-eps(point_val) eps(point_val)]+point_val;
else
    e = [1 2 5 10] * 10^(floor(log10(d)));
    [m i] = min(abs(e - d));
    f = e(i);
    y = f * [floor(min(x)/f - 1) ceil(max(x)/f + 1)];
end
