function [f, v] = poly2fv(x, y)
%POLY2FV Convert polygonal region to patch faces and vertices
%
%   [F, V] = POLY2FV(X, Y) converts the polygonal region represented by the
%   contours (X, Y) into a faces matrix, F, and a vertices matrix, V, that
%   can be used with the PATCH function to display the region. If the
%   polygon represented by X and Y has multiple parts, either the
%   NaN-separated vector format or the cell array format may be used. The
%   POLY2FV function creates triangular faces.
%
%   Most Mapping Toolbox functions adhere to the convention that individual
%   contours with clockwise-ordered vertices are external contours and
%   individual contours with counterclockwise-ordered vertices are internal
%   contours. Although the POLY2FV function ignores vertex order, you
%   should follow the convention when creating contours to ensure
%   consistency with other functions.
%
%   Example
%   -------
%   Display a rectangular region with two holes using a single patch
%   object.
%
%       % External contour, rectangle.
%       x1 = [0 0 6 6 0];
%       y1 = [0 3 3 0 0];
%      
%       % First hole contour, square.
%       x2 = [1 2 2 1 1];
%       y2 = [1 1 2 2 1];
%
%       % Second hole contour, triangle.
%       x3 = [4 5 4 4];
%       y3 = [1 1 2 1];
%
%       % Compute face and vertex matrices.
%       [f, v] = poly2fv({x1, x2, x3}, {y1, y2, y3});
%
%       % Display the patch.
%       patch('Faces', f, 'Vertices', v, 'FaceColor', 'r', ...
%             'EdgeColor', 'none');
%       axis off, axis equal
%
%   See the documentation for POLYBOOL for additional examples illustrating
%   POLY2FV.
%
%   See also ISPOLYCW, PATCH, POLY2CW, POLY2CCW, POLYBOOL.

% Copyright 2004-2010 The MathWorks, Inc.
% $Revision: 1.1.4.4 $  $Date: 2010/03/22 03:51:46 $

if isempty(x)
   f = [];
   v = [];
   return;
end

if iscell(x)
    [x, y] = polyjoin(x,y);
else
    checkxy(x, y, mfilename, 'X', 'Y', 1, 2)
end

p = vectorsToGPC(x, y, mfilename, 'x', 'y');

t = gpcmex('poly2tri', p);

if isempty(t)
   v = [];
   f = [];
   return;
else
   total_num_vertices = numel(t(1).x);
   total_num_faces = total_num_vertices - 2;
   for k = 2:numel(t)
      total_num_vertices = total_num_vertices + numel(t(k).x);
      total_num_faces = total_num_faces + numel(t(k).x) - 2;
   end
end

v = zeros(total_num_vertices, 2);
f = zeros(total_num_faces, 3);

v_offset = 0;
f_offset = 0;
for k = 1:numel(t)
   nv = numel(t(k).x);
   vk = [t(k).x(:) t(k).y(:)];
   fk = (1:(nv - 2)).' + v_offset;
   fk = [fk, fk+1, fk+2];
   v((1:nv) + v_offset, :) = vk;
   f((1:(nv-2)) + f_offset, :) = fk;
   v_offset = v_offset + nv;
   f_offset = f_offset + nv - 2;
end
