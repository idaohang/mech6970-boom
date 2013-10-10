function varargout = behrmann(varargin)
%BEHRMANN  Behrmann Cylindrical Projection
%
%  This is an orthographic projection onto a cylinder secant at the
%  30 degree parallels.  It is equal area, but distortion of shape
%  increases with distance from the standard parallels.  Scale is true
%  along the standard parallels and constant between two parallels
%  equidistant from the Equator.  This projection is not equidistant.
%
%  The projection is named for Walter Behrmann, who presented it in
%  1910 and is a special form of the Equal Area Cylindrical projection
%  secant at 30 degrees N and S.

% Copyright 1996-2008 The MathWorks, Inc.
% $Revision: 1.9.4.5 $  $Date: 2008/06/16 16:47:41 $

if nargin == 1
    % Set defaults.
    mstruct = varargin{1};
    [mstruct.trimlat, mstruct.trimlon, mstruct.mapparallels] ...
        = fromDegrees(mstruct.angleunits, [-90 90], [-180 180], 30);
    mstruct.nparallels = 0;
    mstruct.fixedorient  = [];
    varargout = {mstruct};
else
    varargout = cell(1,max(nargout,1));
    [varargout{:}] = eqacylin(varargin{:});
end
