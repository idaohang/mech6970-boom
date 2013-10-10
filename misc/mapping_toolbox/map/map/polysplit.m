function [latcells,loncells] = polysplit(lat,lon)
%POLYSPLIT Convert line or polygon parts from vector form to cell arrays
%
%   [LATCELLS,LONCELLS] = POLYSPLIT(LAT,LON) returns the NaN-delimited
%   segments of the vectors LAT and LON as N-by-1 cell arrays with one
%   polygon segment per cell.  LAT and LON must be the same size and have
%   identically-placed NaNs.  The polygon segments are column vectors if
%   LAT and LON are column vectors, and row vectors otherwise.
%
%   See also isShapeMultipart, POLYJOIN, POLYBOOL, POLYCUT.

% Copyright 1996-2007 The MathWorks, Inc.
% $Revision: 1.4.4.6 $  $Date: 2007/10/10 20:47:45 $

checkinput(lat,{'numeric'},{'real','vector'},mfilename,'LAT',1);
checkinput(lon,{'numeric'},{'real','vector'},mfilename,'LON',2);
[lat, lon] = removeExtraNanSeparators(lat, lon);

% Find NaN locations.
indx = find(isnan(lat(:)));

% Simulate the trailing NaN if it's missing.
if ~isempty(lat) && ~isnan(lat(end))
    indx(end+1,1) = numel(lat) + 1;
end

%  Extract each segment into pre-allocated N-by-1 cell arrays, where N is
%  the number of polygon segments.  (Add a leading zero to the indx array
%  to make indexing work for the first segment.)
N = numel(indx);
latcells = cell(N,1);
loncells = cell(N,1);
indx = [0; indx];
for k = 1:N
    iStart = indx(k)   + 1;
    iEnd   = indx(k+1) - 1;
    latcells{k} = lat(iStart:iEnd);
    loncells{k} = lon(iStart:iEnd);
end
