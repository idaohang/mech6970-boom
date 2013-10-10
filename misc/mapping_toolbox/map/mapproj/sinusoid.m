function varargout = sinusoid(varargin)
%SINUSOID  Sinusoidal Pseudocylindrical Projection
%
%  This projection is equal area.  Scale is true along every parallel
%  and along the central meridian.  There is no distortion along the
%  Equator or along the central meridian, but it becomes severe near
%  the outer meridians at high latitudes.
%
%  This projection was developed in the 16th century.  It was used
%  by Jean Cossin in 1570 and by Jodocus Hondius in Mercator atlases
%  of the early 17th century.  It is the oldest pseudocylindrical
%  projection currently in use, and is sometimes called the
%  Sanson-Flamsteed or the Mercator Equal Area projections.

% Copyright 1996-2008 The MathWorks, Inc.
% $Revision: 1.8.4.4 $  $Date: 2008/06/16 16:47:56 $

if nargin == 1
    % Set defaults.
    mstruct = varargin{1};
    [mstruct.trimlat, mstruct.trimlon] ...
        = fromDegrees(mstruct.angleunits, [-90 90], [-180 180]);
    mstruct.mapparallels = 0;
    mstruct.nparallels   = 0;
    mstruct.fixedorient  = [];
    varargout = {mstruct};
else
    varargout = cell(1,max(nargout,1));
    [varargout{:}] = bonne(varargin{:});
end
