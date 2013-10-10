function bbox = mapbbox(varargin)
%MAPBBOX Compute bounding box of georeferenced image or data grid
%
%   BBOX = MAPBBOX(R, HEIGHT, WIDTH) computes the 2-by-2 bounding box of
%   a georeferenced image or regular gridded data set in a planar map
%   coordinate system.  R is either a 3-by-2 referencing matrix defining
%   a 2-dimensional affine transformation from intrinsic pixel
%   coordinates to map coordinates, or a spatialref.MapRasterReference
%   object.  HEIGHT and WIDTH are the image dimensions.  BBOX bounds the
%   outer edges of the image in map coordinates:
%
%                           [minX minY
%                            maxX maxY]
%
%   BBOX = MAPBBOX(R, SIZEA) accepts SIZEA = [HEIGHT, WIDTH, ...]
%   instead of HEIGHT and WIDTH.
%
%   BBOX = MAPBBOX(INFO) accepts a scalar struct array with the fields:
%
%              'RefMatrix'   A 3-by-2 referencing matrix
%              'Height'      A scalar number
%              'Width'       A scalar number
%
%   See also GEOTIFFINFO, MAKEREFMAT, MAPOUTLINE, PIXCENTERS, PIX2MAP.

% Copyright 1996-2010 The MathWorks, Inc.
% $Revision: 1.1.6.4 $  $Date: 2010/10/11 14:47:20 $ 

error(nargchk(1,3,nargin,'struct'))

[R, h, w] = parsePixMapInputs('MAPBBOX', [], varargin{:});

if isobject(R)
    % Compute bounding box from spatial referencing object.
    bbox = [sort(R.XLimWorld)' sort(R.YLimWorld)'];
else
    % Compute bounding box from referencing matrix.
    bbox = mapbbox_refmat(R, h, w);
end

%-----------------------------------------------------------------------------

function bbox = mapbbox_refmat(R, h, w)
% Compute bounding box from referencing matrix

outline = [(0.5 + [0  0;...
                   0  w;...
                   h  w;...
                   h  0]), ones(4,1)] * R;

bbox = [min(outline); max(outline)];
