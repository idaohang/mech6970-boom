function [x,y] = closePolygonInRectangle(x, y, xLimit, yLimit, tolSnap)
%closePolygonInRectangle Close polygon parts along edges of rectangle
%
%   [X,Y] = closePolygonInRectangle(X, Y, xLimit, yLimit, tolSnap)
%   connects the end points of open polygon loops that terminate along
%   the edges of the rectangle defined by the two-element vectors xLimit
%   and yLimit. It assumes that xLimit(1) < xLimit(2) and yLimit(1) <
%   yLimit(2). The polygons are represented by the NaN-separated (and
%   terminated) column vectors X and Y. X and Y are assumed to be
%   consistent with a topological convention in which, as one traverses
%   the polygon vertices in the order in which they are stored, polygon
%   interiors fall on the right-hand side of the curve and exteriors (or
%   voids) fall on the left-hand side.  That is, a closed, clockwise
%   ring represents an outer boundary and a closed, counter-clockwise
%   ring encloses an interior hole or void. The scalar tolerance
%   tolSnap, used for snapping the endpoints of open loops onto the
%   perimeter of the rectangle, should be substantially smaller than
%   diff(xLimit) and diff(yLimit).

% Copyright 2008-2010 The MathWorks, Inc.
% $Revision: 1.1.6.4 $  $Date: 2010/06/07 16:34:22 $

if ~isempty(x)
    [x, y] = snapToLimits(x, y, xLimit, yLimit, tolSnap);
    [xTraced, yTraced, isOpen] = traceAllOpenCurves(x, y, xLimit, yLimit);
    [x, y] = removeOpenCurves(x, y, isOpen);
    x = vertcat(xTraced{:}, x);
    y = vertcat(yTraced{:}, y);
end

%--------------------------------------------------------------------------

function [x,y] = snapToLimits(x, y, xLimit, yLimit, tolSnap)
% For each open curve in the NaN-separated arrays (x,y), if the first or
% last vertex is within distance tolSnap of the bounding rectangle defined
% by xLimit and yLimit, snap it onto the boundary.

% Indices of the first and last vertex in each curve.
[first, last] = internal.map.findFirstLastNonNan(x);

% Identify open curves.
isOpen = ~((x(first) == x(last)) & (y(first) == y(last)));

% Snap first vertices that belong to open
% curves and are close to the limits.
x(first(isOpen & (abs(x(first) - xLimit(1)) < tolSnap))) = xLimit(1);
x(first(isOpen & (abs(x(first) - xLimit(2)) < tolSnap))) = xLimit(2);
y(first(isOpen & (abs(y(first) - yLimit(1)) < tolSnap))) = yLimit(1);
y(first(isOpen & (abs(y(first) - yLimit(2)) < tolSnap))) = yLimit(2);

% Snap last vertices that belong to open
% curves and are close to the limits.
x(last(isOpen & (abs(x(last) - xLimit(1)) < tolSnap))) = xLimit(1);
x(last(isOpen & (abs(x(last) - xLimit(2)) < tolSnap))) = xLimit(2);
y(last(isOpen & (abs(y(last) - yLimit(1)) < tolSnap))) = yLimit(1);
y(last(isOpen & (abs(y(last) - yLimit(2)) < tolSnap))) = yLimit(2);

%--------------------------------------------------------------------------

function [xTraced, yTraced, isOpen] ...
    = traceAllOpenCurves(x, y, xlimit, ylimit)
% Link up/close up open curves along the edges of the rectangle,
% returning the results in a pair of cell arrays.

% Compute structures with information on the end points of each open
% curve, including their normalized "positions" as measured clockwise
% from the lower-left corner of the rectangle.
[C, E, isOpen] = organizeOpenCurves(x, y, xlimit, ylimit);

% Loop over each open curve. If it hasn't yet been traced, then trace it
% from start to end, and go on the next open curve (moving clockwise
% around the perimeter of the rectangle -- in order of increasing
% "position") and inserting vertices for any corners traversed along the
% way. Stop when the first vertex of the starting curve is reached.
traced = false(size(C));

% Pre-allocate cell arrays contain to contain the results. (Unused cells
% will be removed later.)
xTraced = cell(size(C));
yTraced = cell(size(C));

% Corner coordinates, indexed in order of increasing position.
xCorner = xlimit([1 1 2 2])';
yCorner = ylimit([1 2 2 1])';

% Trace each open curve exactly once.
n = 0;
for k = 1:numel(C)
    if ~traced(k)
        % Trace starting with the right side (clockwise)
        [xPoly, yPoly, traced] = traceSimplePolygon( ...
            k, traced, x, y, C, E, xCorner, yCorner); 
        n = n + 1;
        xTraced{n} = xPoly;
        yTraced{n} = yPoly;
    end
end

% Remove unused cells.
unused = cellfun(@isempty,xTraced);
xTraced(unused) = [];
yTraced(unused) = [];

%--------------------------------------------------------------------------

function [xPoly, yPoly, traced] = traceSimplePolygon( ...
   k, traced, x, y, C, E, xCorner, yCorner)
% Trace the simple, closed polygon that contains the k-th open curve,
% returning the results in vectors xPoly and yPoly. Update the logical
% array TRACED to keep track of which open curves have been traversed.

% Remove from E all curves that have already been
% traced and used to form other polygons.
E(traced([E.Index])) = [];

% Remove last points from E, keeping only start points.
E(~[E.IsFirst]) = [];

% Save copy of E to reuse on second trace.
Esaved = E;

% Number of untraced curves
nCurves = numel(E);

% Trace once to count the vertices.
numVertices = 0;
j = k;
done = false;
n = 0;
while ~done
    % Ensure that the loop cannot run forever.
    assert(n < nCurves, ...
        'map:closePolygonInRectangle:tracingFailed1', ...
        'Tracing failed to converge in first loop.')
    n = n + 1;
    
    numVertices = numVertices + (C(j).Last - C(j).First + 1);
    endPosition = C(j).PosLast;
    
    [j, corners, E] = findNextCurve(endPosition, E);
    
    numVertices = numVertices + numel(corners);
    done = (j == k);
end

% Two extra elements: Repeat first vertex followed by terminating NaN. 
numVertices = numVertices + 2;

% Trace again to assign the vertices.
E = Esaved;
xPoly = zeros(numVertices,1);
yPoly = zeros(numVertices,1);
j = k;
i1 = 1;
done = false;
n = 0;
while ~done
    % Ensure that the loop cannot run forever.
    assert(n < nCurves, ...
        'map:closePolygonInRectangle:tracingFailed2', ...
        'Tracing failed to converge in second loop.')
    n = n + 1;

    endPosition = C(j).PosLast;        
    subs = (C(j).First):(C(j).Last);
    traced(j) = true;
    i2 = i1 + C(j).Last - C(j).First;
    xPoly(i1:i2) = x(subs);
    yPoly(i1:i2) = y(subs);
    
    [j, corners, E] = findNextCurve(endPosition, E);
    
    i1 = i2 + 1;
    i2 = i2 + numel(corners);
    xPoly(i1:i2) = xCorner(corners);
    yPoly(i1:i2) = yCorner(corners);
    i1 = i1 + numel(corners);
    done = (j == k);
end

if (numVertices > 2) ...
        && ((xPoly(end-2) ~= xPoly(1)) || (yPoly(end-2) ~= yPoly(1)))
    % Repeat the first vertex, if necessary.
    xPoly(end-1) = xPoly(1);
    yPoly(end-1) = yPoly(1);
else
    % Otherwise, since the pre-allocation of xPoly and yPoly provided
    % slots in which to replicate xPoly(1) and yPoly(1), these arrays
    % were over sized by one element each and need to be shortened
    % accordingly.
    xPoly(end) = [];
    yPoly(end) = [];
end

% Append NaN-terminators.
xPoly(end) = NaN;
yPoly(end) = NaN;

%--------------------------------------------------------------------------

function [j, corners, E] = findNextCurve(endPosition, E)
% Traversing the edge of the rectangle in a clockwise direction from a
% point with "position" equal to endPosition, find the closest start
% point and return the corresponding "C" structure index, j. In
% addition, return a list of the  corners, if any, that must be
% traversed when following the edge of the rectangle clockwise to the
% start point of curve j. Update E to remove the start point so that it
% can't be found the next time around.

p = find(endPosition <= [E.Position]);
if isempty(p)
    p = 1;
else
    p = p(1);
end
j = E(p).Index;
corners = cornersTraversed(endPosition, E(p).Position);
E(p) = [];

%--------------------------------------------------------------------------

function corners = cornersTraversed(position1, position2)
% Return, in sequence, the indices of the corners that are traversed when
% moving from position 1 to position 2 in clockwise order (order of
% increasing position).

% All possible corners in increasing order
corners = (0:7);  
    
% Unwrap positions into the half-open interval [0 8).
if position2 < position1
    position2 = position2 + 4;
end

% Eliminate the corners that are not traversed.
corners(corners <= position1) = [];
corners(corners >= position2) = [];

% Wrap back to into the interval [0 3], and add one to convert from
% position to index.  The final result contains up to four contiguous
% elements from the array [1 2 3 4 1 2 3].
corners = 1 + mod(corners,4);

%--------------------------------------------------------------------------

function [x, y] = removeOpenCurves(x, y, remove)
% Given NaN-separated vertex arrays x and y and a logical vector REMOVE
% which has one element per curve, remove the curves for which REMOVE is
% true.

[first, last] = internal.map.findFirstLastNonNan(x);
first(~remove) = [];
last( ~remove) = [];
for k = 1:numel(first)
    x(first(k):last(k)) = NaN;
    y(first(k):last(k)) = NaN;
end
[x, y] = removeExtraNanSeparators(x, y);
