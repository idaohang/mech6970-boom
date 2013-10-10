function checkxy(x, y, function_name, x_var_name, y_var_name, x_pos, y_pos)
%CHECKXY Check validity of map x and y vectors
%
%   CHECKXY(X, Y, FUNCTION_NAME, X_VAR_NAME, Y_VAR_NAME, X_POS, X_POS)
%   ensures that X and Y are real vectors of matching size and equal NaN
%   locations.

% Copyright 2006 The MathWorks, Inc.
% $Revision: 1.1.6.2 $  $Date: 2006/06/15 20:10:55 $

% Input arguments are not checked for validity.

% Check numeric, 2d, and real.
checkinput(x, ...
   {'numeric'}, {'real','2d','vector'}, function_name, x_var_name, x_pos);
checkinput(y, ...
   {'numeric'}, {'real','2d','vector'}, function_name, y_var_name, y_pos);

if ~isequal(isnan(x), isnan(y))
   eid = sprintf('%s:%s:inconsistentXY', getcomp, function_name);  
   msg = sprintf(...
      'Function %s expected its %s and %s input arguments,\n%s and %s, %s', ... 
      upper(function_name), num2ordinal(x_pos), num2ordinal(y_pos), ...
      x_var_name, y_var_name, 'to match in size or NaN locations.');
   error(eid, '%s', msg);
end
