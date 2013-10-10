function [xg, yg] = gluePolygonsOnVerticalEdges(x, y, xLimit, tolSnap)
% Wrap an infinite vertical strip in the plane into a cylinder, gluing
% together polygons that meet on opposing vertical edges. x and y are
% vertex arrays for topologically valid polygons, stored as
% NaN-separated column vectors. Another way to think about this: Take a
% set of polygons bounded by a rectangle and glue the left and right
% edges together to form a cylinder. Except for NaN-separators,
%
%    all(xLimit(1) <= x <= xLimit(2)).
%
% Assume input polygons are closed to within tolSnap.

% Copyright 2010 The MathWorks, Inc.
% $Revision: 1.1.6.2 $  $Date: 2010/06/07 16:34:09 $

if all(isnan(x))
    x = [];
    y = [];
end

if isempty(x)
    xg = [];
    yg = [];
    return
end

[x, y] = removeExtraNanSeparators(x, y);

% Snap values of x that fall almost on the line x == xLimit(1) or
% x == xLimit(2).
[x, y] = snapToLimits(x, y, xLimit, tolSnap);

% Unless the polygons extend to both limits, there's no work to do.
if (xLimit(1) < min(x)) || (max(x) < xLimit(2))
    xg = x;
    yg = y;
    return
end

% Allocate xg and yg to hold the "glued" polygons. As a result of
% gluing, the number of vertices should decrease, so the original size of
% x and y should be more than enough. There are two advantages to
% starting with NaN-filled vectors: (1) We can drop in parts without
% having to explicitly add NaN-separators, we just skip a slot between
% subsequent parts and (2) we can clear out unused slots when done by
% calling removeExtraNanSeparators (although that does not turn out to
% be necessary here.)
xg = NaN(size(x));
yg = NaN(size(x));

% If any vertical strip along the edge x == xLimit(1) overlaps with a
% vertical strip along the edge x == xLimit(2), remove the overlapping
% portion from both edges.
[x, y, first, last] = preprocessEdges(x, y, xLimit);

% At this point, every curve that touches one (or both) of the vertical
% edges and that needs to be traced from one limit to the other is open,
% because it no longer has a segment along that edge to close it.

% Identify open curves.
isOpen = ~((x(first) == x(last)) & (y(first) == y(last)));
firstOpen = first(isOpen);
lastOpen  = last(isOpen);

% Each open curve should end at (or very close to) the start point of
% another open curve, if we consider the points (xlimit(1), y) and
% (xlimit(2), y) to be identical. Construct a lookup table which, given
% an open curve index k returns the index of the open curve next(k)
% whose start point coincides with curve k's end point.
next = nextCurveLookup(x, y, xLimit, tolSnap, firstOpen, lastOpen);

% Trace the open curves, copying vertices into xg and yg.  n, a
% positive, scalar integer indicates where to start when copying
% addition vertices (from closed curves) into xg and yg.
[xg, yg, n] = traceOpenCurves(x, y, firstOpen, lastOpen, next, xg, yg);

% Copy the closed curves.
firstClosed = first(~isOpen);
lastClosed  = last(~isOpen);
[xg, yg, n] = copyClosedCurves(x, y, firstClosed, lastClosed, xg, yg, n);

% Removed excess NaNs from the ends of xg and yg.
xg(n:end) = [];
yg(n:end) = [];

%-----------------------------------------------------------------------

function [x, y, first, last] = preprocessEdges(x, y, xLimit)

% Remove extra vertices from the edges x == xLimit(1) and x == xLimit(2).
% A vertex is "extra" if both the preceding and following vertices are
% on the same edge.
[x, y] = removeExtraEdgeVertices(x, y, xLimit(1));
[x, y] = removeExtraEdgeVertices(x, y, xLimit(2));

% Identify segments that run along the left and right edges (and remove
% any redundant vertices that may be present within those segments).
% The outputs leftPair and rightPairs are 2-by-n arrays in indices in
% which each column corresponds to a single pair.
[leftPairs,  x, y] = polygonEdgePairs(x, y, xLimit(1));
[rightPairs, x, y] = polygonEdgePairs(x, y, xLimit(2));

% The coordinate arrays below, yLeft and yRight, will also be 2-by-n.
yLeft  = y(leftPairs);
yRight = y(rightPairs);

[first, last] = internal.map.findFirstLastNonNan(x);

% Remove from each segment the part(s) that overlap with segments on the
% opposite edge. Do this by identifying a break point (the index from
% the top row) and a set of points (with NaN-separators and terminators)
% to insert. Loop over the columns in leftPairs ...
for k = 1:size(leftPairs,2)
    ySeg = nonOverlap(yLeft(:,k), yRight([2 1],:));
    [x, y, first, last] = removeOverlap( ...
        x, y, first, last, leftPairs(:,k), ySeg, xLimit(1));
end

% Repeat on the right.
for k = 1:size(rightPairs,2)
    ySeg = -nonOverlap(-yRight(:,k), -yLeft([2 1],:));
    [x, y, first, last] = removeOverlap( ...
        x, y, first, last, rightPairs(:,k), ySeg, xLimit(2));
end

% Re-sort indices of start and end points; remove degenerate vertices.
first = sort(first);
last  = sort(last);
degenerate = (first == last);
first(degenerate) = [];
last( degenerate) = [];

%-----------------------------------------------------------------------

function [x, y, first, last] = removeOverlap( ...
    x, y, first, last, pair, ySeg, xBound)
% Given the indices to a pair of vertices on a vertical edge (PAIR, a
% 2-element column vector) and set of intervals in y (ySeg), remove from
% the multi-part lines x and y the sections that between the vertices in
% the PAIR that correspond to intervals in which the interval defined by
% sorting y(pair) overlaps with ySeg.

noOverlap = isequal(ySeg(1:end-1), y(pair));
if ~noOverlap
    if isempty(ySeg)
        % Designate a break by adding indices to first and last. Do
        % not modify the vertex arrays themselves.
        last( end+1,1) = pair(1);
        first(end+1,1) = pair(2);
    else
        % Designate a break and append to the vertex arrays x and y one
        % or more segments that run along the edge without fully
        % covering the interval y(pair).
        [x, y, first, last] = appendSegments( ...
            x, y, first, last, pair, ySeg, xBound);        
    end
end

%-----------------------------------------------------------------------

function [x, y, first, last] ...
    = appendSegments(x, y, first, last, edgePair, ySeg, xBound)
% Append segments that connect pairs of vertices along the vertical
% edges of a rectangle.

first(end+1,1) = edgePair(2);
last( end+1,1) = edgePair(1);

n = find(isnan(ySeg));
segFirst = numel(x) + [1, 1 + n(1:n-1)];
segLast  = numel(x) + n - 1;

first = [first segFirst];
last =  [last  segLast];

xSeg = NaN(size(ySeg));
xSeg(~isnan(ySeg)) = xBound;

x = [x; xSeg];
y = [y; ySeg];

%--------------------------------------------------------------------------

function [x,y] = snapToLimits(x, y, xLimit, tolSnap)
% For each open curve in the NaN-separated arrays (x,y), if the first or
% last vertex is within distance tolSnap of the limits defined
% by xLimit, snap it to the boundary.

% Adapted from subfunction in private/closePolygonInRectangle.

% Indices of the first and last vertex in each curve.
[first, last] = internal.map.findFirstLastNonNan(x);

% Identify open curves.
isOpen = ~((x(first) == x(last)) & (y(first) == y(last)));

% Snap first vertices that belong to open
% curves and are close to the limits.
x(first(isOpen & (abs(x(first) - xLimit(1)) < tolSnap))) = xLimit(1);
x(first(isOpen & (abs(x(first) - xLimit(2)) < tolSnap))) = xLimit(2);

% Snap last vertices that belong to open
% curves and are close to the limits.
x(last(isOpen & (abs(x(last) - xLimit(1)) < tolSnap))) = xLimit(1);
x(last(isOpen & (abs(x(last) - xLimit(2)) < tolSnap))) = xLimit(2);

%--------------------------------------------------------------------------

function next = nextCurveLookup(x, y, xLimit, tolSnap, firstOpen, lastOpen)
% Array for looking up index of curve starting where k-th one ends.

% Start points as column vectors
xs = x(firstOpen);
ys = y(firstOpen);

% End points as row vectors
xe = x(lastOpen)';
ye = y(lastOpen)';

% Expand into n-by-n arrays.
n = numel(firstOpen);
nOnes = ones(1,n);

xs = xs(:,nOnes);
ys = ys(:,nOnes);

xe = xe(nOnes, :);
ye = ye(nOnes, :);

% n-by-n logical connects is true when an end point coincides
% (or nearly coincides) with a start point.
connects = abs(ye - ys) <= tolSnap & (abs(xe - xs) <= tolSnap ...
    | (xs == xLimit(1) & xe == xLimit(2)) ...
    | (xs == xLimit(2) & xe == xLimit(1)));

[next, ~] = find(connects);

%--------------------------------------------------------------------------

function [xg, yg, n] ...
    = traceOpenCurves(x, y, firstOpen, lastOpen, next, xg, yg)
% Trace the open curves, copying vertices into xg and yg. From each end
% point, look for a (nearly) coincident start point and go from there.
% Use a logical array to keep track of which curves have already been
% traced.  Return an index value n that is equals 2 plus the index of
% the last vertex copied into xg and yg. In other words, n is the place
% to start when copying addition vertices (from closed curves) into xg
% and yg.

% Keep track of the number of curves traced in relation to the total
% number, to ensure that the while loop terminates (see assert below).
nOpen = numel(firstOpen);
nTraced = 0;

n = 1;
traced = false(size(firstOpen));
for k = 1:numel(firstOpen)
    if ~traced(k)
        % Start tracing a set of open curves that form a closed ring on
        % the sphere.
        
        % Copy vertex 1 from the k-th open curve.
        xg(n) = x(firstOpen(k));
        yg(n) = y(firstOpen(k));
        
        % Copy the rest of the k-th open curves and the curves to which
        % it connects.
        j = k;
        while ~traced(j);
            % Ensure that the loop cannot run forever. The following
            % assertion should never be triggered.
            nTraced = nTraced + 1;
            assert(nTraced <= nOpen, ...
                'map:gluePolygonsOnVerticalEdges:tracingFailed', ...
                'Failed to converge when tracing open curves.')
            
            % Copy vertices 2 through end from the j-th open curve.
            n = n + 1;
            s = 1 + firstOpen(j);
            e = lastOpen(j);
            m = n + e - s;
            xg(n:m) = x(s:e);
            yg(n:m) = y(s:e);
            
            % Set up for next iteration
            n = m;
            traced(j) = true;
            j = next(j);
        end
        
        % Advance by 2 instead of 1, leaving a NaN-separator in xg and yg to
        % separate this curve from the next one.
        n = n + 2;        
    end
end

%--------------------------------------------------------------------------

function [xg, yg, n] ...
    = copyClosedCurves(x, y, firstClosed, lastClosed, xg, yg, n)
% Copy the closed curves, taking the index-to-next-vertex-position, n,
% as an input and returning an updated value.

for k = 1:numel(firstClosed)
    s = firstClosed(k);
    e = lastClosed(k);
    % Compute m such that n:m has the same size as s:e
    m = n + e - s;
    xg(n:m) = x(s:e);
    yg(n:m) = y(s:e);
    
    % Advance by 2 instead of 1, leaving a NaN-separator in xg and yg to
    % separate this curve from the next one.
    n = m + 2;
end
