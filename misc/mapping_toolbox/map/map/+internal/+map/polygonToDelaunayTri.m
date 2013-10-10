function dt = polygonToDelaunayTri(x, y)
% Construct a Delaunay Triangulation object for a polygon represented by
% vertex arrays X and Y.  The polygon may comprise multiple,
% NaN-separated parts, including inner rings that outline voids.  Outer
% rings have their vertices ordered in a clockwise direction and inner
% rings are counter-clockwise.  (For both inner and outer rings, the
% inside of the polygon is always on the right-hand side when traversing
% the vertices in the order given.)

% Copyright 2009-2010 The MathWorks, Inc.
% $Revision: 1.1.6.1 $  $Date: 2010/05/13 17:37:04 $

% Work with column vectors.
x = x(:);
y = y(:);

% Repair and/or simplify the polygon topology, if necessary.
[x, y] = polybool('union', x, y, x, y);

% Locate NaN-separators in the vertex arrays (x and y have the
% same size and have NaNs in the corresponding locations).
n = isnan(x(:));

% Remove NaN-separators from the vertex arrays.
x(n) = [];
y(n) = [];

% Build an edge-constraints matrix for use with DelaunayTri.
constraints = buildConstraintsMatrix(n);

% Construct a constrained Delaunay triangulation object. Turn off
% warnings about adjustment of vertices.
w(3) = warning('off','MATLAB:DelaunayTri:ConsSplitPtWarnId');
w(2) = warning('off','MATLAB:DelaunayTri:DupPtsConsUpdatedWarnId');
w(1) = warning('off','MATLAB:DelaunayTri:ConsConsSplitWarnId');
c = onCleanup(@() warning(w));
dt = DelaunayTri(x, y, constraints);

%-----------------------------------------------------------------------

function c = buildConstraintsMatrix(n)

% Given N, a logical array having the same size as the vertex arrays and
% containing true for each NaN-valued element of the vertex arrays and
% false for all others, build C, a (number of edges)-by-2 edge
% constraints matrix for use with DelaunayTri.  C contains indices which
% refer to elements of the vertex arrays in the positions they have once
% all the NaN-separators are removed.  It has one row per edge.  The
% first element in the row is an index for the start point of that edge,
% and the second element in the row is an index for the endpoint of that
% edge.

% Number of elements in vertex arrays, including NaN-separators
m = numel(n);

% Ignore for the moment any NaN-separators that may be present and
% construct a column vector of start-point indices.  It begins like this:
% i1 = [1 2 3 4 5 6 ...]'
i1 = (1:(m-1))';

% Continuing to ignore Nan-separators, construct a column of end-point
% indices.  It begins like this:  i2 = [2 3 4 5 6 7 ...]'
i2 = (2:m)';

% Combine i1 and i2 to construct a preliminary version of the
% edge-constraints matrix.  It begins like this:
%
%    c = [1  2;
%         2  3;
%         3  4;
%         4  5;
%         5  6;
%         6  7;
%         ...]
%
% Note that at this point, if there are no NaN-separators (there's just
% one simple ring, that is), then we're done.  We've just linked each
% vertex to the next one traversing our way around the ring.
c = [i1 i2];

% Next, correct c to account for the presence of NaN-separators, if any.
% Here's the first correction:  Without changing its size, reduce the
% values in c to account for the fact that the Nan-separators will be
% removed from the vertex arrays before they are passed to DelaunayTri.
% Before reaching the first separator, no corrections are required.
% After the reaching the first, reduce the values in c by 1. Then, after
% reaching second, reduce the values in c by 2. And so on.  
%
% If there are NaN-separators in positions 5 and 10, for example, then
% we'd have:
%                n = [0 0 0 0 1 0 0 0 0 1 0 0 ...]',
%
%    t = cumsum(n) = [0 0 0 0 1 1 1 1 1 2 2 2 ...]', and
%
%    c = [1  2;
%         2  3;
%         3  4;
%         4  4;
%         4  5;
%         5  6;
%         6  7;
%         7  8;
%         8  8;
%         8  9;
%         9 10;
%         ...]
% 
t = cumsum(double(n));
c = c - [t(i1) t(i2)];

% At this point c contains indices that are valid for the vertex
% arrays with their NaN-separators removed, but it contains extraneous
% rows that attempt to link the various rings across the separators.
% In the example, the 4-th row [4 4] and 5-th row [4 5] are both
% extraneous. [4 4] attempts to connect the last vertex in the first ring
% to itself, and [4 5] attempts to connect the last vertex in the first
% ring to the first vertex in the second ring. Instead, there should be
% a clean break between the end of one ring and the start of the next,
% like this:
%
%    c = [1  2;
%         2  3;
%         3  4;
%         5  6;
%         6  7;
%         7  8;
%         9 10;
%         ...]
%
%  Thus the second and final correction is the elimination of any
%  extraneous rows. As the example suggests, these rows come in pairs.
%  The index of the first row in each pair corresponds to the indices of
%  the NaNs and the index of the second row exceeds the index of the
%  first by 1, so we want just need to eliminate each row for which the
%  logical column vector n(i1) | n(i2) contains true.
c(n(i1) | n(i2), :) = [];
