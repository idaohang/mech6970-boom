function [lat,lon] = interpm(lat,lon,maxdiff,method,units)
%INTERPM  Densify latitude-longitude sampling in lines or polygons
%
%  [lat,lon] = INTERPM(lat,long,maxdiff) linearly interpolates between
%  vector data coordinate points where necessary to return data with no
%  two connected points separated by an angular distance greater than
%  maxdiff. Maxdiff must be in the same units as the input lat and lon
%  data.
%
%  [lat,lon] = INTERPM(lat,long,maxdiff,'method') interpolates between
%  vector data coordinate points using a specified interpolation method.
%  Valid interpolation methods strings are 'gc' for great circle, 'rh'
%  for rhumb lines, and 'lin' for linear interpolation. With no units
%  specified, lat,long and maxdiff are assumed to be in units of degrees.
%
%  [lat,lon] = INTERPM(lat,long,maxdiff,'method','units') interpolates
%  between vector data coordinate points using a specified interpolation
%  method. Inputs and outputs are in the specified units.
%
%  See also INTRPLAT, INTRPLON, RESIZEM.

% Copyright 1996-2007 The MathWorks, Inc.
% $Revision: 1.10.4.4 $  $Date: 2007/11/09 20:24:12 $

error(nargchk(3, 5, nargin, 'struct'))

if nargin == 5
    if ~ischar(method)
        error(['map:' mfilename ':mapError'], 'Method must be a string')
    end
    if ~ischar(units)
        error(['map:' mfilename ':mapError'], 'Units must be a string')
    end
elseif nargin == 4
    if ~ischar(method)
        error(['map:' mfilename ':mapError'], 'Method must be a string')
    end
    units = 'deg';
elseif nargin == 3;
    method = 'lin';
    units = 'deg';
end

%  Dimensional tests

if ~isequal(size(lat),size(lon))
    error(['map:' mfilename ':mapError'], ...
        'Inconsistent lat and lon dimensions')
elseif max(size(maxdiff)) ~= 1
    error(['map:' mfilename ':mapError'], ...
        'Scalar maximum angular difference required')
elseif any([~isreal(lat) ~isreal(lon) ~isreal(maxdiff)])
    wid = sprintf('%s:%s:ignoreComplex', getcomp, mfilename);
    warning(wid, 'Imaginary parts of complex arguments ignored')
	lat = real(lat);
    lon = real(lon);
    maxdiff = real(maxdiff);
end

%  Test the track string

if isempty(method)
    method = 'lin';       %  Default is linear interpolation
else
    validstr = {'gc','rh','lin'};
    indx = find(strcmpi(method,validstr));
    if numel(indx) ~= 1
        error(['map:' mfilename ':mapError'], 'Unrecognized track string')
    else
        method = validstr{indx};
    end
end

[lat, lon] = doInterpm(lat,lon,maxdiff,method,units);
