function [value, varargs] = findNameValuePair(name, default, varargin)
%FINDNAMEVALUEPAIR Parse name-value pair and return value
%
%   [VALUE, VARARGS] = findNameValuePair(NAME, DEFAULT, VARARGIN) returns
%   the value of the name-value pair. NAME is a string. DEFAULT is the
%   value of VALUE if NAME is not found. VARARGIN is the input to search
%   for the name string. VALUE is the argument in VARARGIN after NAME.
%   VARARGS is VARARGIN with the name/value pair removed if found.
%
%   See also PARSEPV.

% Copyright 2009 The MathWorks, Inc.
% $Revision: 1.1.6.1 $  $Date: 2009/05/14 17:05:40 $

% Set value to default for not-found case.
value = default;

% Set the output varargs to the input varargins.
varargs = varargin;

% Check varargin for name, value pairs, if found, delete them and set the
% value.
if ~isempty(varargin)
   deleteIndex = false(size(varargin));
   internal.map.checkNameValuePairs(varargin{:});

   for k = 1:2:numel(varargin)
      try
         validatestring(varargin{k}, {name}, mfilename);
         value = varargin{k+1};
         deleteIndex(k:k+1) = true;
      catch e
         % If validatestring threw an error, we can ignore it because all
         % it means is that we didn't find the name we were looking for,
         % but shouldn't ignore unexpected errors.
         firstPart = ['MATLAB:' mfilename ':'];
         if ~strncmp(e.identifier, firstPart, numel(firstPart))
             % Error is not from validatestring.
             rethrow(e)
         end
      end
   end
   varargs(deleteIndex)=[];
end
