function numArgs = getNumberOfDataArgs(varargin)
%GETNUMBEROFDATAARGS Return number to first string input
%
%   NUMARGS = getNumberOfDataArgs(VARARGIN) returns the number of arguments
%   preceding the first string-valued argument. If VARARGIN is empty or no
%   string arguments are found, NUMARGS is 0.
%
%   See also PARSEPV.

% Copyright 2009 The MathWorks, Inc.
% $Revision: 1.1.6.1 $  $Date: 2009/05/14 17:05:41 $

numArgs = nargin;
for i=1:nargin
   if ischar(varargin{i})
      numArgs = i-1;
      return;
   end
end
