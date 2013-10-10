function [C, E, isOpenOnBoundary] = organizeOpenCurves(x, y, xLimit, yLimit)
% Given a set of curves defined by the NaN-separated vectors X and Y and
% bounded by the rectangle defined by XLIMIT and YLIMIT, compute several
% properties of the subset of curves which are open. (X,Y) may comprise a
% combination of simple closed curves and open curves, with the added
% constraint that every open curve begins and ends on one of the edges of
% the rectangle.
%
% The key property computed here is the "position" of each end point.
% Position is defined to be a normalized distance from the lower left
% corner (where x == xLimit(1) and y == yLimit(1)), measured clockwise
% along the perimeter of the rectangle. Along the horizontal edges,
% distance is normalized by the width of the rectangle, diff(xLimit) and
% along the vertical edges it is normalized by the height of the rectangle,
% diff(yLimit).
%
% Because of the normalization scheme, the corner points of the rectangle
% have integer position values:
%
%    Lower left corner:   Position = 0
%    Upper left corner:   Position = 1
%    Upper right corner:  Position = 2
%    Lower right corner:  Position = 3
%
% The position property thus makes it possible to identify corners and to
% tell if one must pass one or more corners when traveling along the edge
% of the rectangle from one of the endpoints of one of the open curves to
% one of the endpoints of a different curve (or the other endpoint of the
% same curve).
%
% This function returns two structure arrays. The first output, C, contains
% one element for each of the open curves and has four fields:
%
%    C(k).First -- The index of the first vertex of the k-th open curve in
%                  the coordinate arrays x and y
%
%    C(k).Last  -- The index of the last vertex of the k-th open curve in
%                  the coordinate arrays x and y
%
%    C(k).PosFirst -- Position value of the first vertex of the k-th open
%                     curve
%
%    C(k).PosLast -- Position value of the last vertex of the k-th open
%                    curve
%
% The second output, E, has one element for each open-curve end point. (So
% it has twice as many elements as C.) E has the three fields:
%
%    E(k).Index -- Indexes the corresponding element of C
%
%    E(k).IsFirst -- True for the start point of an open curve, false for
%                    the end point
%
%    E(k).Position -- Position value of the corresponding end point
%
% Array E is sorted in terms of ascending values of position.
%
% Also returned isOpenOnBoundary, a logical array with one element for
% each ring or open curve defined by the NaN-separated vertex arrays x
% and y. For each part, the corresponding element of isOpenOnBoundary is
% true if and only if its first vertex identical differs from its
% last vertex and both fall on the boundary of the rectangle.

% Copyright 2008-2009 The MathWorks, Inc.
% $Revision: 1.1.6.2 $  $Date: 2009/10/16 05:00:45 $

% Indices of the first and last vertex in each curve.
[first, last] = internal.map.findFirstLastNonNan(x);

% Identify open curves that terminate on the boundary of the rectangle.
isOpen = ~((x(first) == x(last)) & (y(first) == y(last)));
isOpenOnBoundary = isOpen & ...
    onBoundary(x(first), y(first), xLimit, yLimit) & ...
    onBoundary(x(last),  y(last),  xLimit, yLimit);
openCurves = find(isOpenOnBoundary);

if ~any(isOpenOnBoundary)
    C = [];
    E = [];
else
    % Initialize C, with one element per open curve.
    C(numel(openCurves),1) = struct( ...
        'First', [], ...
        'Last', [], ...
        'PosFirst', [], ...
        'PosLast', []);
    
    % Store first and last vertex indices in C.
    for k = 1:numel(openCurves)
        j = openCurves(k);
        C(k).First = first(j);
        C(k).Last  = last(j);
    end
    
    % Initialize E, a structure array of endpoints.
    E(2*numel(openCurves),1) = struct( ...
        'Index', [], ...
        'IsFirst', [], ...
        'Position', []);
    
    % Allocate endpoint coordinate arrays to use in computing position.
    ex = zeros(size(E));
    ey = zeros(size(E));
    
    % Assign Index and IsFirst fields and ex and ey arrays.
    for k = 1:2:(2*numel(openCurves))
        index = 1 + (k-1)/2;
        E(k).Index = index;
        E(k).IsFirst = true;
        E(k+1).Index = index;
        E(k+1).IsFirst = false;
        j = openCurves(index);
        ex(k) = x(first(j));
        ey(k) = y(first(j));
        ex(k+1) = x(last(j));
        ey(k+1) = y(last(j));
    end
    
    % Assign position values in E, processing the four edges in sequence.
    position = zeros(numel(E),1);
    q = (ey == yLimit(1));
    position(q) = 3 + (xLimit(2) - ex(q))/diff(xLimit);
    q = (ex == xLimit(2));
    position(q) = 2 + (yLimit(2) - ey(q))/diff(yLimit);
    q = (ey == yLimit(2));
    position(q) = 1 + (ex(q) - xLimit(1))/diff(xLimit);
    q = (ex == xLimit(1));
    position(q) = 0 + (ey(q) - yLimit(1))/diff(yLimit);
    for k = 1:numel(E)
        E(k).Position = position(k);
    end
    
    % Assign position values in C (for both end points of each curve).
    for k = 1:numel(C)
        C(k).PosFirst = position(2*k-1);
        C(k).PosLast  = position(2*k);
    end
    
    % Sort E by position.
    [~, indx] = sort(position);
    E = E(indx);
end

%--------------------------------------------------------------------------

function tf = onBoundary(x, y, xLimit, yLimit)
% True if and only if (x,y) falls on the boundary of the rectangle
% defined by limits xLimit and yLimit.

% tf = ...
%     (x == xLimit(1)) | ...
%     (x == xLimit(2)) | ...
%     (y == yLimit(1)) | ...
%     (y == yLimit(2));

tf = ...
    abs(x - xLimit(1)) < eps(10) | ...
    abs(x - xLimit(2)) < eps(10) | ...
    abs(y - yLimit(1)) < eps(10) | ...
    abs(y - yLimit(2)) < eps(10);
