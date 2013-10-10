function [lat,lon,alt] = minvtran(varargin)
%MINVTRAN Unproject features from map to geographic coordinates
%
%   [LAT, LON] = MINVTRAN(X, Y) applies the inverse transformation
%   defined by the map projection in the current map axes, converting
%   point locations and line/polygon vertices in a planar, projected map
%   coordinate system to latitudes and longitudes.
%
%   [LAT, LON, HEIGHT] = MINVTRAN(X, Y, Z) applies the inverse
%   projection to 3-D input, resulting in 3-D output.  If the input
%   Z is empty or omitted, then Z = 0 is assumed.
%
%   [...] = MINVTRAN(MSTRUCT,...) takes a valid map projection structure
%   as the first argument. In this case, no map axes is needed.
%
%   See also MAPS, MFWDTRAN, PROJINV, PROJLIST.

% Copyright 1996-2008 The MathWorks, Inc.
% $Revision: 1.1.6.6 $  $Date: 2008/06/16 16:47:54 $

% Additional syntax for internal use only:
%
%   [LAT, LON, HEIGHT] = MINVTRAN(X, Y, Z, objectType, SAVEPTS)
%   removes all clips and trims from the input data. The input SAVEPTS
%   is a structure, returned by function MFWDTRAN, which contains
%   information regarding how the object was clipped and trimmed by a
%   previous, forward transformation. Allowable objectType strings are:
%
%       'surface' for map graticules
%       'line' for line objects
%       'patch' for patch objects
%       'light' for light objects
%       'text' for text objects
%       'none' to ignore any clipping and trimming of input data.

[mstruct, x, y, z, objectType, savepts] = parseInputs(varargin{:});
[lat, lon, alt] = feval(mstruct.mapprojection, ...
    mstruct, x, y, z, objectType, 'inverse', savepts);

%-----------------------------------------------------------------------

function [mstruct, x, y, z, objectType, savepts] = parseInputs(varargin)

if (nargin >= 1) && isstruct(varargin{1})
    error(nargchk(3, 6, nargin, 'struct'))
    mstruct = varargin{1};
    varargin(1) = [];
else
    error(nargchk(2, 5, nargin, 'struct'))
    mstruct = gcm;
end

if mstruct.geoid(1) <= 0;
    error(['map:' mfilename ':mapprojError'], ...
    'Semimajor axis of reference ellipsoid must be nonzero and positive.')
end

x = varargin{1};
y = varargin{2};

% Assign default values as needed.
defaults = { ...
    zeros(size(x)), ...                  % z
    'none', ...                          % objectType
    struct('trimmed',[],'clipped',[])};  % savepts

varargin(end+1:numel(defaults)+2) = defaults(numel(varargin)-1:end);

z          = varargin{3};
objectType = varargin{4};
savepts    = varargin{5};

% Ensure non-empty z even in the case where varargin{3} is []
if isempty(z)
    z = zeros(size(x));
end
