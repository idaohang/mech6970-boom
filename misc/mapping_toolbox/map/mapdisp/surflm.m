function h = surflm(lat, lon, Z, varargin)
%SURFLM 3-D shaded surface with lighting on map axes
%
%   SURFLM(LAT, LON, Z) and SURFLM(LATLIM, LONLIM, Z) are the same as
%   SURFM(...) except that they highlight the surface with a light
%   source. The default light source (45 degrees counterclockwise from
%   the current view) and reflectance constants are the same as in
%   SURFL.
%
%   SURFLM(...,S) and SURFLM(...,S,K) use a light source vector, S, and a
%   vector of reflectance constants, K.  See the help for SURFL for more
%   information on S and K.
%
%   H = SURFLM(...) returns a handle to the surface object.
%
%   See also SURFACE, SURFACEM, SURFL, SURFM.

% Copyright 1996-2009 The MathWorks, Inc.
% $Revision: 1.15.4.7 $  $Date: 2009/05/14 17:06:29 $

if nargin == 0
    % Error on obsolete syntax SURFLM
    errorOnMissingUI(mfilename)
elseif (nargin == 1) || (nargin == 2)
    % Error on obsolete syntaxes SURFLM(Z) and SURFLM(Z,GRATSIZE)
    errorOnObsoleteSyntax(nargin)
end

if numel(lat) == 2 && numel(lon) == 2
    % SURFLM(LATLIM, LONLIM, Z, ...)
    latlim = lat;
    lonlim = lon;
    [lat,lon] = meshgrat(latlim,lonlim,size(Z));
else
    % SURFLM(LAT, LON, Z, ...)
    if ~isequal(size(lat),size(lon),size(Z))
        error(['map:' mfilename ':sizeMismatch'], ...
            'LAT, LON, and Z must have consistent sizes.')
    end
end

%  Validate map axes
ax = getParentAxesFromArgList(varargin);
if ~isempty(ax)
    mstruct = gcm(ax);
else
    mstruct = gcm;
end

%  Project the surface data
if ~strcmp(mstruct.mapprojection,'globe')
    [x,y,ignored,savepts] = mfwdtran(mstruct,lat,lon,[],'surface');
else
    error(['map:' mfilename ':axesIsGlobe'], ...
        'SURFLM cannot be used with a GLOBE map axes.');
end

%  Display the map
nextmap(varargin);
h0 = surfl(x,y,Z,varargin{:});
otherprops = {...
    'ButtonDownFcn', 'uimaptbx',...
    'UserData',       savepts,...
    'EdgeColor',     'none'};
set(h0,otherprops{:})

%  Set handle return argument if necessary
if nargout > 0;
    h = h0;
end

%-----------------------------------------------------------------------

function errorOnObsoleteSyntax(argcount)

% Issue an error in case someone tries one of the following obsolete
% syntaxes:
%
% SURFLM(Z) and SURFLM(Z,S) project the data grid by constructing
% a map graticule to span the MapLatLimit and MapLonLimit specified
% in the map structure.

if argcount == 1
    syntax = 'Z';
    eid = ['map:' mfilename ':obsoleteSyntax1'];
else
    syntax = 'Z,S';
    eid = ['map:' mfilename ':obsoleteSyntax2'];
end

error( eid, ...
    ['The syntax', ...
    '\n\n   surflm(%s)\n\n', ...
    'made the implicit assumption of a perfect match between the', ...
    ' geographic\nlimits of the current map axes', ...
    ' and your data grid, and this syntax is\nno longer supported.' ...
    ' Instead, you need to pass the geographic limits\nof your data', ...
    ' to SURFLM, like this:', ...
    '\n\n   surflm(latlim,lonlim,%s)\n', ...
    '\nwhere LATLIM and LONLIM have the form:\n\n', ...
    '   LATLIM = [southern_limit northern_limit]\n', ...
    '   LONLIM = [western_limit  eastern_limit ]'], ...
    syntax, syntax)
