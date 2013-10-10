function checkrefvec(refvec, func_name, var_name, arg_pos)
%CHECKREFVEC Check validity of referencing vector
%
%     This function is intentionally undocumented and is intended for
%     use only by other Mapping Toolbox functions.  Its behavior may
%     change, or the function itself may be removed, in a future
%     release.
%
%   CHECKREFVEC(REFVEC, FUNC_NAME, VAR_NAME, ARG_POS) ensures that the
%   referencing vector REFVEC is a 1-by-3 vector of real-valued doubles
%   and that the first element is positive.
%
%   See also CHECKREFMAT.

% Copyright 2007-2009 The MathWorks, Inc.
% $Revision: 1.1.6.3 $  $Date: 2009/03/09 19:16:13 $

validateattributes(...
    refvec, {'double'}, {'real','finite','nonempty','vector'}, ...
    func_name, var_name, arg_pos);

% REFVEC must have three elements.
assert(numel(refvec) == 3, ...
    ['map:' func_name ':refvecNumelNot3'], ...
    ['Function %s expected its %s input argument, %s,\n' ...
    'to have 3 elements.'], ...
    upper(func_name), num2ordinal(arg_pos), var_name)

% The first element of REFVEC (cells/angleunit) must be positive.
assert(refvec(1) > 0, ...
    ['map:' func_name ':cellsPerAngleUnitNotPositive'], ...
    'Function %s expected the first element of %s\nto be positive.', ...
    upper(func_name), var_name)
