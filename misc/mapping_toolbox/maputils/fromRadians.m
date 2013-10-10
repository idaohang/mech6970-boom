function varargout = fromRadians(toUnits, varargin)
%fromRadians Convert angles from radians
%
%   [angle1, angle2, ...]
%       = fromRadians(toUnits, angle1InRadians, angle2InRadians, ...)
%   converts angle1InRadians, angle2InRadians, ... from radians to the
%   specified output ("to") angle units.  toUnits can be either
%   'degrees' or 'radians' and may be abbreviated.  The inputs
%   angle1InRadians, angle2InRadians, ... and their corresponding
%   outputs are numeric arrays of various sizes, with size(angleN)
%   matching size(angleNInRadians).
%
%   See also: fromDegrees, radtodeg, toDegrees, toRadians.

% Copyright 2009 The MathWorks, Inc.
% $Revision: 1.1.6.1 $  $Date: 2009/03/09 19:31:50 $

varargout = abstractAngleConv( ...
    'radians', 'degrees', @radtodeg, toUnits, varargin{:});
