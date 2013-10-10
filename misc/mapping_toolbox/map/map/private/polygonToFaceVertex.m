function [faces, vertices] = polygonToFaceVertex(x,y)
%polygonToFaceVertex  Face-vertex decomposition of multipart polygon
%
%   Uses an implementation based on constrained Delaunay triangulation.
%
%   See also POLY2FV

% Copyright 2009-2010 The MathWorks, Inc.
% $Revision: 1.1.6.2 $  $Date: 2010/05/13 17:37:16 $

% Construct a constrained Delaunay triangulation.
dt = internal.map.polygonToDelaunayTri(x, y);

% Get the vertices and faces.
vertices = dt.X;
faces = dt.Triangulation;

% Discard faces that fall outside the polygon.
inside = dt.inOutStatus();
faces(~inside,:) = [];
