function h = pcolorm(lat, lon, Z, varargin)
%PCOLORM Project regular data grid on map axes in z = 0 plane
%
%   PCOLORM(LAT, LON, Z) constructs a surface to represent the data grid
%   Z in the current map axes.  The surface lies flat in the horizontal
%   plane with its 'CData' property set to Z. LAT and LON are vectors or
%   2-D arrays that define the latitude-longitude graticule mesh on
%   which Z is displayed. See SURFM for a complete description of the
%   various forms that LAT and LON can take. PCOLORM will clear the
%   current map if the hold state is 'off'.
%
%   PCOLORM(LATLIM, LONLIM, Z) defines the graticule using the latitude
%   and longitude limits LATLIM and LONLIM, which should match the
%   geographic extent of the data grid Z.  LATLIM is a two-element
%   vector of the form:
%
%               [southern_limit northern_limit]
% 
%   Likewise, LONLIM has the form:
%
%                [western_limit eastern_limit]
%
%   A latitude-longitude graticule of size 50-by-100 is constructed. The
%   surface 'FaceColor' property is 'texturemap' except when Z is
%   precisely 50-by-100, in which case it is 'flat'.
%
%   PCOLORM(..., PROP1, VAL1, PROP2, VAL2,...) applies additional
%   MATLAB graphics properties to the surface, via property/value pairs.
%   Any property accepted by the SURFACE function may be specified,
%   except for XData, YData, and ZData.
%
%   H = PCOLORM(...) returns a handle to the surface object.
%
%   See also GEOSHOW, MESHM, SURFACEM, SURFM.

% Copyright 1996-2008 The MathWorks, Inc.
% $Revision: 1.9.4.7 $  $Date: 2008/05/12 21:33:08 $

if nargin == 1
    % Error on obsolete syntax PCOLORM(Z)
    errorOnObsoleteSyntax
elseif nargin == 2
    % Error on obsolete syntax PCOLORM(Z,GRATSIZE)
    errorOnObsoleteSurfmSyntax(mfilename,2)
end
error(nargchk(3, Inf, nargin, 'struct'))

% The approach here is to always call SURFACEM with a syntax that will
% result in ALT being an array of zeros.  This happens if we restrict
% ourselves to the following syntaxes:
%
% H = SURFACEM(LAT, LON, Z)
% H = SURFACEM(LAT, LON, Z, <property/value pairs>)
nextmap;
if nargout == 0
    surfacem(lat, lon, Z, varargin{:}); 
else
    h = surfacem(lat, lon, Z, varargin{:});
end

%-----------------------------------------------------------------------

function errorOnObsoleteSyntax

% Dedicated function for the obsolete syntax PCOLORM(Z). 

error(['map:' mfilename ':obsoleteSyntax1'], ...
    ['The syntax', ...
    '\n\n   pcolorm(Z)\n\n', ...
    'made the implicit assumption of a perfect match between the', ...
    ' geographic\nlimits of the current map axes', ...
    ' and your data grid, and this syntax is\nno longer supported.' ...
    ' Instead, you need to call MESHGRAT before calling\nPCOLORM,', ...
    ' like this:', ...
    '\n\n   [lat, lon] = meshgrat(latlim, lonlim, size(Z));', ...
    '\n   pcolorm(lat, lon, Z)\n', ...
    '\nwhere LATLIM and LONLIM have the form:\n\n', ...
    '   LATLIM = [southern_limit northern_limit]\n', ...
    '   LONLIM = [western_limit  eastern_limit ]'])
