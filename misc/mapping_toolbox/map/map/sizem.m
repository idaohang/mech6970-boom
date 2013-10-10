function [r,c,refvec] = sizem(latlim,lonlim,scale)
%SIZEM  Row and column dimensions needed for regular data grid
%
%   [R,C] = SIZEM(LATLIM,LONLIM,SCALE) computes the row and column
%   dimensions needed for a regular data grid aligned with geographic
%   coordinates.  LATLIM and LONLIM are two-element vectors defining the
%   latitude and longitude limits in degrees. SCALE is a scalar specifying
%   the number of data samples per unit of latitude and longitude (e.g. 10
%   entries per degree).
%
%   SZ = SIZEM(...) returns a single output, where SZ = [R C].
%
%   [R,C,REFVEC] = SIZEM(...) returns the referencing vector for the data
%   grid.

% Copyright 1996-2007 The MathWorks, Inc.
% $Revision: 1.9.4.5 $  $Date: 2007/11/09 20:25:32 $
% Written by:  E. Byrns, E. Brown

%  Argument tests

if ~isequal(sort(size(latlim)),sort(size(lonlim)),[1 2])
    error(['map:' mfilename ':mapError'], ...
        'Lat and long limits must be 2 element vectors')
elseif max(size(scale)) ~= 1
    error(['map:' mfilename ':mapError'], ...
        'Scale input must be scalar')
end

latlim  = ignoreComplex(latlim,  mfilename, 'latlim');
lonlim  = ignoreComplex(lonlim,  mfilename, 'lonlim');
scale   = ignoreComplex(scale,   mfilename, 'scale');

%  Determine the starting and ending latitude and longitude

startlat = min(latlim);
endlat   = max(latlim);
startlon = lonlim(1);
endlon   = lonlim(2);
if endlon < startlon
    endlon = endlon + 360;
end

%  Compute the number of rows and columns needed

rows = ceil((endlat - startlat)*scale);
cols = ceil((endlon - startlon)*scale);

%  Set the output arguments

if nargout == 1
    r = [rows cols];
elseif nargout == 2
    r = rows;   c = cols;
elseif nargout == 3
    r = rows;   c = cols;  refvec = [scale endlat startlon];
end
