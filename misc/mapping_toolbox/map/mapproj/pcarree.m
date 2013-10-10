function varargout = pcarree(varargin)
%PCARREE  Plate Carree Cylindrical Projection
%
%  This is a projection onto a cylinder tangent at the Equator.
%  Distortion of both shape and area increase with distance from
%  the Equator.  Scale is true along all meridians (i.e., it is
%  equidistant) and the Equator, and is constant along any parallel
%  and along the parallel of opposite sign.
%
%  This projection, like the more general equidistant cylindrical, is
%  credited to Marinus of Tyre, thought to have invented it about
%  A.D. 100.  It may, in fact, have been originated by Eratosthenes,
%  who lived approximately 275-195 B.C.  The Plate Carree has the
%  most simply constructed graticule of any projection.  It was
%  used frequently in the 15th and 16th centuries, and is quite common
%  today in very simple computer mapping programs.  It is the simplest
%  and limiting form of the equidistant cylindrical projection.  Another
%  name for this projection is the Simple Cylindrical.  Its transverse
%  aspect is the Cassini projection.

% Copyright 1996-2008 The MathWorks, Inc.
% $Revision: 1.10.4.4 $  $Date: 2008/06/16 16:47:55 $

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
    [varargout{:}] = eqdcylin(varargin{:});
end
