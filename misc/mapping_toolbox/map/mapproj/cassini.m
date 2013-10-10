function varargout = cassini(varargin)
%CASSINI  Cassini Transverse Cylindrical Projection
%
%  This is a projection onto a cylinder tangent at the central
%  meridian.  Distortion of both shape and area are functions of distance
%  from the central meridian.  Scale is true along the central meridian
%  and along any straight line perpendicular to the central meridian
%  (i.e., it is equidistant).
%
%  This projection is the transverse aspect of the Plate Carree projection,
%  developed by Cesar Francois Cassini de Thury (1714-84).  It is still
%  used for the topographic mapping of a few countries.
%
%  See also CASSINISTD.

% Copyright 1996-2008 The MathWorks, Inc.
% $Revision: 1.8.4.5 $  $Date: 2008/06/16 16:47:43 $

if nargin == 1
    % Set defaults.
    mstruct = varargin{1};
    [mstruct.trimlat, mstruct.trimlon, mstruct.fixedorient] ...
        = fromDegrees(mstruct.angleunits, [-90 90], [-180 180], -90);
    mstruct.mapparallels = 0;
    mstruct.nparallels   = 0;
    varargout = {mstruct};
else
    varargout = cell(1,max(nargout,1));
    [varargout{:}] = eqdcylin(varargin{:});
end
