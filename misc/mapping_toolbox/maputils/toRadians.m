function varargout = toRadians(fromUnits, varargin)
%toRadians Convert angles to radians
%
%   [angle1InRadians, angle2InRadians, ...]
%       = toRadians(fromUnits, angle1, angle2, ...)
%   converts angle1, angle2, ... to radians from the specified input
%   ("from") angle units.  fromUnits can be either 'degrees' or 'radians'
%   and may be abbreviated.  The inputs angle1, angle2, ... and their
%   corresponding outputs are numeric arrays of various sizes, with
%   size(angleNInRadians) matching size(angleN).
%
%   See also: degtorad, fromDegrees, fromRadians, toDegrees.

% Copyright 2009 The MathWorks, Inc.
% $Revision: 1.1.6.1 $  $Date: 2009/03/09 19:31:54 $

varargout = abstractAngleConv( ...
    'radians', 'degrees', @degtorad, fromUnits, varargin{:});
