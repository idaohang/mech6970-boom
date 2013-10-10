function tf = multipleDistinctLocations(lat,lon,angleunits)
% Return true if and only if the lat, lon arrays include more than one
% distinct location.

% Copyright 2007 The MathWorks, Inc.
% $Revision: 1.1.6.1 $  $Date: 2007/05/10 13:48:04 $

% Merge lat, lon into a single complex array and see if it contains more
% than one unique, finite value.
lon = mod(toDegrees(angleunits, lon), 360);
s = lat(:) + i*lon(:);
s(~isfinite(s)) = [];
tf = (numel(unique(s)) > 1);
