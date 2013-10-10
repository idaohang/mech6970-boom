function [tf, wrapsOrigin] = ispolycwPolar(th, r)
%ispolycw Polar coordinate counterpart to ispolycw
%
%   TF is true if vertices are in clockwise order, for a multipart
%   polygon in polar coordinates. Vertex vectors TH and R has the same
%   meaning as in CART2POL and POL2CART. wrapsOrigin is true if the
%   polygon wraps the origin (in either a clockwise or counterclockwise
%   direction).
%
%   Note that one _cannot_ simply apply pol2cart and then work in a
%   Cartesian system. The reason is that a linear segment between two
%   vertices that is a straight line in range-azimuth will map to a curve
%   in the Cartesian system, so to do that we'd have to first interpolate
%   extra vertices and select a threshold for a sufficiently dense
%   sampling, etc., and all results would still be approximate.
%
%   See also ISPOLYCW.

% Copyright 2009 The MathWorks, Inc.
% $Revision: 1.1.6.2 $  $Date: 2010/01/19 02:55:27 $

[first, last] = internal.map.findFirstLastNonNan(th);
th = unwrapMultipart(th);
tf = false(size(first));
wrapsOrigin = tf;
for k = 1:numel(first)
    thk = th(first(k):last(k));
    rk = r(first(k):last(k));
    accumulatedAngle = sum(diff(thk));
    % We expect accumulatedAngle to come out very close to one of these
    % three values:
    %
    %  2*pi ==> Counter-clockwise curve that wraps once around the origin
    %     0 ==> Closed curve that does not wrap
    % -2*pi ==> Clockwise curve that wraps once around the origin
    if accumulatedAngle > pi
        % In the interval (pi Inf]; result must be close to 2*pi.
        tf(k) = false;
        wrapsOrigin(k) = true;
    elseif accumulatedAngle < -pi
        % In the interval [-Inf -pi); result must be close to -2*pi.
        tf(k) = true;
        wrapsOrigin(k) = true;
    else
        % In the interval [-pi pi]; result must be close to 0.
        tf(k) = ~ispolycw(thk, rk);
        wrapsOrigin(k) = false;
    end
end
