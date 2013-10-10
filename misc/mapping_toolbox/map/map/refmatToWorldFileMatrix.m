function W = refmatToWorldFileMatrix(refmat)
%refmatToWorldFileMatrix Convert referencing matrix to world file matrix
%
%   W = refmatToWorldFileMatrix(REFMAT) converts the 3-by-2 referencing
%   matrix REFMAT to a 2-by-3 world file matrix W.
%
%   For the definition of a referencing matrix, see the help for
%   MAKEREFMAT.
%
%   For the definition of a world file matrix for a planar system,
%   see the help for the worldFileMatrix method of the
%   spatialref.MapRasterReference class.
%
%   For the definition of a world file matrix for a geographic system,
%   see the help for the worldFileMatrix method of the
%   spatialref.GeoRasterReference class.
%
%   See also MAKEREFMAT, worldFileMatrixToRefmat
%      spatialref.GeoRasterReference/worldFileMatrix
%      spatialref.MapRasterReference/worldFileMatrix

% Copyright 2010 The MathWorks, Inc.
% $Revision: 1.1.6.1.2.1 $  $Date: 2010/12/03 21:43:41 $

% The following expressions are derived in worldFileMatrixToRefmat.m

internal.map.validateRasterReference(refmat, ...
    {}, 'refmatToWorldFileMatrix', 'REFMAT', 1)

Cinv = [0  1  1;...
        1  0  1;...
        0  0  1];

W = refmat' * Cinv;
