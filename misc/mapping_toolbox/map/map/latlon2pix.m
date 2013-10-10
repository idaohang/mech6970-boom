function [row, col] = latlon2pix(R, lat, lon)
%LATLON2PIX Convert latitude-longitude coordinates to pixel coordinates
%
%   [ROW, COL] = LATLON2PIX(R,LAT,LON) calculates pixel  coordinates ROW,
%   COL from latitude-longitude coordinates LAT, LON.  R is either a 3-by-2
%   referencing matrix that transforms intrinsic pixel coordinates to
%   geographic coordinates, or a spatialref.GeoRasterReference object.  LAT
%   and LON are vectors or arrays of matching size.  The outputs ROW and
%   COL have the same size as LAT and LON.  LAT and LON must be in degrees.
%
%   Longitude wrapping is handled: Results are invariant under the
%   substitution LON = LON +/- N * 360 where N is an integer.  Any point on
%   the earth that is included in the image or gridded data set
%   corresponding to R will yield row/column values between 0.5 and 0.5 +
%   the image height/width, regardless of what longitude convention is
%   used.
%
%   Example 
%   -------
%      % Find the pixel coordinates of the upper left and lower right 
%      % outer corners of a 2-by-2 degree gridded data set.  
%      R = makerefmat( 1, 89, 2, 2);
%      [UL_row, UL_col] = latlon2pix(R,  90, 0)     % Upper left
%      [LR_row, LR_col] = latlon2pix(R, -90, 360)   % Lower right
%      [LL_row, LL_col] = latlon2pix(R, -90, 0)     % Lower left
%      % Note that the in both the 2nd case and 3rd case we get a column
%      % value of 0.5, because the left and right edges are on the same
%      % meridian and (-90, 360) is the same point as (-90, 0).
%
%   See also MAKEREFMAT, PIX2LATLON, MAP2PIX.

% Copyright 1996-2010 The MathWorks, Inc.
% $Revision: 1.1.10.8 $  $Date: 2010/10/11 14:47:15 $ 

error(nargchk(3,3,nargin,'struct'))

% Validate referencing matrix or GeoRasterReference object.
internal.map.validateRasterReference(R, ...
    {'spatialref.GeoRasterReference'}, 'latlon2pix', 'R', 1)

if isobject(R)
    [col, row] = R.geographicToIntrinsic(lat, lon);
else
    [row, col] = latlon2pix_refmat(R, lat, lon);
end

%--------------------------------------------------------------------------

function [row, col] = latlon2pix_refmat(R, lat, lon)
% Compute the transformation using referencing matrix R

cycle = 360;  % Degrees only for now

% Resolve longitude ambiguity:
% For which values of n, if any is row >= 0.5 and col >= 0.5, where
%   [row, col] = map2pix(R, lon + cycle * n, lat)?

% Start with values for n = 0 and n = 1 (all we need because of linearity).
[row0, col0] = map2pix(R, lon, lat);          % n = 0
[row1, col1] = map2pix(R, lon + cycle, lat);  % n = 1

% Find limiting values of n as separately constrained by the rows and
% columns.
[rLower, rUpper] = findLimits(row0,row1);
[cLower, cUpper] = findLimits(col0,col1);

% Choose a value for n within the intersection of the limits (if possible)
n = max(rLower,cLower);
t = min(rUpper,cUpper);
n(n == -Inf) = t(n == -Inf);
n(n ==  Inf) = 0;

[row, col] = map2pix(R, lon + cycle * n, lat);

%--------------------------------------------------------------------------

function [lowerLim, upperLim] = findLimits(c0, c1)

d = c1 - c0;
Z = (0.5 - c0) ./ d;

lowerLim = -Inf * ones(size(Z));
lowerLim(d > 0) = ceil(Z(d > 0));

upperLim = Inf * ones(size(Z));
upperLim(d < 0) = floor(Z(d < 0));
