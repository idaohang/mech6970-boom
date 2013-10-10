function refmat = convertToReferencingMatrix(R)
% Convert GeoRasterReference object to referencing matrix

% Copyright 2009-2010 The MathWorks, Inc.
% $Revision: 1.1.6.2 $  $Date: 2010/09/24 14:33:47 $

deltaLat = R.DeltaLatNumerator / R.DeltaLatDenominator;
deltaLon = R.DeltaLonNumerator / R.DeltaLonDenominator;

refmat = [ ...
    0          deltaLon   R.FirstCornerLon - deltaLon / 2; ...
    deltaLat   0          R.FirstCornerLat - deltaLat / 2]';
