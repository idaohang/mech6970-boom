function refmat = worldFileMatrixToRefmat(W)
%worldFileMatrixToRefmat Convert world file matrix to referencing matrix
%
%   REFMAT = worldFileMatrixToRefmat(W) converts the 2-by-3 world file
%   matrix W to a 3-by-2 referencing matrix REFMAT.
%
%   For the definition of a world file matrix for a planar system,
%   see the help for the worldFileMatrix method of the
%   spatialref.MapRasterReference class.
%
%   For the definition of a world file matrix for a geographic system,
%   see the help for the worldFileMatrix method of the
%   spatialref.GeoRasterReference class.
%
%   For the definition of a referencing matrix, see the help for
%   MAKEREFMAT.
%
%   See also MAKEREFMAT, refmatToWorldFileMatrix,
%      spatialref.GeoRasterReference/worldFileMatrix
%      spatialref.MapRasterReference/worldFileMatrix

% Copyright 2010 The MathWorks, Inc.
% $Revision: 1.1.6.1.2.1 $  $Date: 2010/12/03 21:43:47 $

% World File Matrix to Referencing Matrix Conversion
% --------------------------------------------------
% An affine transformation that is expressed like this in terms of a world
% file matrix:
%
%               [xw yw]' = W * [(xi - 1) (yi - 1) 1]' 
%
% is expressed as:
%
%                     [xw yw] = [yi xi 1] * R
%
% in terms of a referencing matrix R, To obtain R from W, note that
%
%                   [xi-1 yi-1 1]' = C * [yi xi 1]',
%
% where
% 
%                         C = [0  1  -1
%                              1  0  -1
%                              0  0   1].
% 
% Therefore [xw yw]' = W * C * [yi xi 1]'.  Transposing both sides gives
% 
%                     [xw yw] = [yi xi 1] * R
% with
%                           R = (W * C)'.

% Referencing Matrix to World File Matrix Conversion
% --------------------------------------------------
% To reverse the conversion and create a world file from R, use
% 
%                          W = R' * inv(C)
% 
% and
% 
%                     inv(C) = [0  1  1
%                               1  0  1
%                               0  0  1].

W = validateWorldFileMatrix(W, 'worldFileMatrixToRefmat', 'W', 1);

C = [0  1  -1;...
     1  0  -1;...
     0  0   1];

refmat = (W * C)';
