function varargs = designateAxesArgAsParentArg(fcnName, varargin)
%designateAxesArgAsParentArg Designate axes argument as Parent argument
%
%   VARAGS = DesignateAxesArgAsParentArg(fcnName, VARARGIN) returns
%   the VARAGIN cell array with the first argument moved to the end of 
%   VARARGIN preceded by 'Parent' if the first argument is an axes handle.
%   fcnName is just a string containing the function name to be used only
%   in creating error messages.
%
%   See also GEOSHOW, MAPSHOW.

% Copyright 2006-2010 The MathWorks, Inc.
% $Revision: 1.1.6.4 $  $Date: 2010/09/24 14:33:48 $

varargs = varargin;

handleMayBePresent = ...
   ~isempty(varargin) && ...
   isscalar(varargin{1}) && ...
   ishghandle(varargin{1}) && ...
   ~ishghandle(varargin{1},'root');

numDataArgs = internal.map.getNumberOfDataArgs(varargin{:});
firstArgIsHgHandle = handleMayBePresent && ...
  (numDataArgs == 1 || ...                          % GEOSHOW(H,'FILENAME')
   numDataArgs == 2 && isstruct(varargin{2}) || ... % GEOSHOW(H,S)
   numDataArgs > 2);                                % GEOSHOW(H,LAT,LON)
 
if firstArgIsHgHandle
   if ishghandle(varargin{1},'axes')
      if nargin > 2
         % Designate first argument (axes handle) as Parent.
         varargs = [varargin(2:end), {'Parent'}, varargin(1)];
      else
         % Remove the axes from varargin.
         % Let the function show with no data.
         varargs = {};
      end
   else
       error(['map:' fcnName ':handleNotAxesType'], ...
           'Function %s expected its first input to be an ''%s'' handle.', ...
           'axes', upper(fcnName))
   end
end
