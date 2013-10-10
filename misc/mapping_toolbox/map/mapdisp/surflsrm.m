function [hout, msg] = surflsrm(lat,long,map,s,rgbs,clim)
%SURFLSRM 3-D lighted shaded relief of geolocated data grid
%
%  SURFLSRM(lat,long,map) displays the general matrix map colored according
%  to elevation and surface slopes.  By default, shading is based on a light
%  to the east (90 deg.)  at an elevation of 45 degrees.  Also by default,
%  the colormap is constructed from 16 colors and 16 grays.  Lighting is
%  applied before the data is projected.  The current axes must have a valid
%  map projection definition.
%
%  SURFLSRM(lat,long,map,[azim elev]) displays the general matrix map with
%  the light coming from the specified azimuth and elevation.  Angles are
%  specified in degrees, with the azimuth measured clockwise from North,
%  and elevation up from the zero plane of the surface.
%
%  SURFLSRM(lat,long,map,[azim elev],cmap) displays the general matrix map
%  using the provided colormap.  The number of grayscales is chosen to keep
%  the size of the shaded colormap below 256. If the vector of azimuth and
%  elevation is empty, the default locations are used. Color axis limits are
%  computed from the data.
%
%  SURFLSRM(lat,long,map,[azim elev],cmap,clim) uses the provided caxis limits.
%
%  h = SURFLSRM(...) returns the handle to the surface drawn.
%
%  See also MESHLSRM, SHADEREL, MESHM, SURFLM, SURFM, SURFACEM, PCOLORM.

% Copyright 1996-2007 The MathWorks, Inc.
% $Revision: 1.12.4.4 $  $Date: 2007/11/09 20:29:14 $
% Written by:  A. Kim, W. Stumpf

% Obsolete syntax
% ---------------
% [h,msg] = SURFLSRM(...) returns a string indicating any error encountered
if nargout > 1
    warnObsoleteMSGSyntax(mfilename)
    msg = '';
end

%  Initialize outputs
if nargout ~= 0;
    hout = [];
end

%  Input argument tests
if nargin == 0;
    errorOnMissingUI(mfilename)
elseif nargin==3;
    rgbs = [];
    clim = [];
    s = [];
elseif nargin==4;
    rgbs = [];
    clim = [];
elseif nargin==5;
    clim = [];
elseif nargin<=2;
    error(['map:' mfilename ':mapdispError'], ...
        'Incorrect number of arguments')
end

%  Input dimension tests
if any([ndims(lat) ndims(long) ndims(map)] > 2)
        error(['map:' mfilename ':mapdispError'], ...
            'Input lat, long and map matrices must not contain any pages');
elseif ~isequal(size(lat),size(long),size(map))
        error(['map:' mfilename ':mapdispError'], ...
            'Inconsistent dimensions on input lat, lon and map matrices');
end

%  Get the current map structure
mstruct = gcm;

%  Set the light source azimuth and elevation
if ~isempty(s) && length(s) ~= 2
        error(['map:' mfilename ':mapdispError'], ...
            'Light source vector must consist of azimuth and elevation');
end

%  Set the color axis limits
if isempty(clim)
    clim = [min(map(:))   max(map(:))];
elseif length(clim) ~= 2
        error(['map:' mfilename ':mapdispError'], ...
            'Color limits must be a two element vector');
end

%  Build shaded relief colormap
if isempty(rgbs);
    [rgbs,clim] = demcmap(map);
end

[rgbindx,rgbmap,clim] = shaderel(long,lat,map,rgbs,s,[],clim);

%  Display shaded relief map
h = surfacem(lat,long,rgbindx,map); colormap(rgbmap)
caxis(clim)

%  Set handle return argument if necessary
if nargout ~= 0;
    hout = h;
end
