function checkNameValuePairs(varargin)
%CHECKNAMEVALUEPAIRS Check name-value pairs 
%
%   checkNameValuePairs(VARARGIN) checks and validates that VARARGIN 
%   consists of name-value pairs. If not, an error is issued.

% Copyright 2009 The MathWorks, Inc.
% $Revision: 1.1.6.2 $  $Date: 2009/11/09 16:25:42 $

if ~isempty(varargin)
   if rem(length(varargin),2)
      error('map:checkNameValuePairs:invalidPairs', ...
          'The parameter/value inputs must be pairs.')
   end
   
   params = varargin(1:2:end);
   for k=1:length(params)
      if ~ischar(params{k})
         error('map:checkNameValuePairs:invalidParameterString', ...
            'The %s parameter/value pair must be a string followed by value.', ...
            num2ordinal(k))
      end
   end
end
