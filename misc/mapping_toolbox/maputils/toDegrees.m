function varargout = toDegrees(fromUnits, varargin)
%toDegrees Convert angles to degrees
%
%   [angle1InDegrees, angle2InDegrees, ...]
%       = toDegrees(fromUnits, angle1, angle2, ...)
%   converts angle1, angle2, ... to degrees from the specified input
%   ("from") angle units.  fromUnits can be either 'degrees' or 'radians'
%   and may be abbreviated.  The inputs angle1, angle2, ... and their
%   corresponding outputs are numeric arrays of various sizes, with
%   size(angleNInDegrees) matching size(angleN).
%
%   See also: fromDegrees, fromRadians, radtodeg, toRadians.

% Copyright 2009 The MathWorks, Inc.
% $Revision: 1.1.6.1 $  $Date: 2009/03/09 19:31:53 $

varargout = abstractAngleConv( ...
    'degrees', 'radians', @radtodeg, fromUnits, varargin{:});
