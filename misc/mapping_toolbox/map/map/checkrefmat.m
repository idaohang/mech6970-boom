function checkrefmat(refmat, func_name, var_name, arg_pos)
%CHECKREFMAT Check validity of referencing matrix
%
%     This function is intentionally undocumented and is intended for
%     use only by other Mapping Toolbox functions.  Its behavior may
%     change, or the function itself may be removed, in a future
%     release.
%
%   CHECKREFMAT(REFMAT, FUNC_NAME, VAR_NAME, ARG_POS) ensures that the
%   referencing matrix REFMAT is a 3-by-2 matrix of real-valued finite
%   doubles.
%
%   See also CHECKREFVEC.

% Copyright 2007-2009 The MathWorks, Inc.
% $Revision: 1.1.6.1 $  $Date: 2009/03/09 19:16:12 $

validateattributes( ...
    refmat, {'double'} ,{'real','2d','finite','nonempty'}, ...
    func_name, var_name, arg_pos)

% REFMAT must be 3-by-2.
assert(isequal(size(refmat),[3,2]), ...
    ['map:' func_name ':refmatNot3by2'], ...
    ['Function %s expected its %s input argument, %s,\n', ...
    'to have size [3,2].'], ...
    upper(func_name), num2ordinal(arg_pos), var_name);
