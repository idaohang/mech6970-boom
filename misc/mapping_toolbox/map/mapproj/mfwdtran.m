function [x,y,z,savepts] = mfwdtran(varargin)
%MFWDTRAN  Project geographic features to map coordinates
%
%   [X, Y] = MFWDTRAN(LAT, LON) applies the forward transformation
%   defined by the map projection in the current map axes, converting point
%   locations and line/polygon vertices given in latitude and longitude
%   to a planar, projected map coordinate system.
%
%   [X, Y, Z] = MFWDTRAN(LAT, LON, HEIGHT) applies the forward projection
%   to 3-D input, resulting in 3-D output.  If the input HEIGHT is empty or
%   omitted, then HEIGHT = 0 is assumed.
%
%   [...] = MFWDTRAN(MSTRUCT,...) takes a valid map projection structure
%   as the first argument. In this case, no map axes is needed.
%
%  See also MAPS, MINVTRAN, PROJFWD, PROJLIST.

% Copyright 1996-2008 The MathWorks, Inc.
% $Revision: 1.1.6.7 $  $Date: 2008/07/09 18:11:40 $

% Additional syntax for internal use only:
%
%   [X, Y, Z, SAVEPTS] = MFWDTRAN(..., objectType) clips and trims the
%   object during the projection process.  The output SAVEPTS is a
%   structure which is associated with each object displayed. It
%   contains information regarding the clips and trims associated with
%   the object. Allowable objectType strings are:
%
%       'surface' for map graticules
%       'line' for line objects
%       'patch' for patch objects
%       'light' for light objects
%       'text' for text objects
%       'none' to omit any clipping and trimming of input data.

[mstruct, lat, lon, alt, objectType] = parseInputs(varargin{:});
[x, y, z, savepts] = feval( ...
    mstruct.mapprojection, mstruct, lat, lon, alt, objectType, 'forward');

%-----------------------------------------------------------------------

function [mstruct, lat, lon, alt, objectType] = parseInputs(varargin)

if (nargin >= 1) && isstruct(varargin{1})
    error(nargchk(3, 5, nargin, 'struct'))
    mstruct = varargin{1};
    varargin(1) = [];
else
    error(nargchk(2, 4, nargin, 'struct'))
    mstruct = gcm;
end

if mstruct.geoid(1) <= 0
    error(['map:' mfilename ':mapprojError'], ...
        'Semimajor axis of reference ellipsoid must be nonzero and positive.')
end

lat = varargin{1};
lon = varargin{2};

% Assign default values as needed.
defaults = { ...
    zeros(size(lat)), ...  % alt
    'none'};               % objectType

varargin(end+1:numel(defaults)+2) = defaults(numel(varargin)-1:end);

alt        = varargin{3};
objectType = varargin{4};

% Ensure non-empty alt even in the case where varargin{3} is []
if isempty(alt)
    alt = zeros(size(lat));
end
