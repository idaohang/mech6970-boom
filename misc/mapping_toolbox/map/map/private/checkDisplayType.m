function varargs = checkDisplayType(geometry, fcnName, varargin)
%CHECKDISPLAYTYPE Check the DisplayType parameter
%
%   VARARGS = CHECKDISPLAYTYPE(GEOMETRY, FCNNAME, VARARGIN) finds and
%   checks the 'DisplayType' parameter with a given GEOMETRY. If the
%   GEOMETRY does not match the DisplayType parameter, a warning is issued.
%   VARARGS contains VARARGIN with the 'DisplayType' parameter and value
%   removed.

% Copyright 2006-2009 The MathWorks, Inc.
% $Revision: 1.1.6.2 $  $Date: 2009/11/09 16:25:53 $

% Find the DisplayType value.
[displayGeometry, varargs] = ...
   internal.map.findNameValuePair('DisplayType', geometry, varargin{:});

% If the displayGeometry from the command line does not match a required
% geometry value, issue a warning.
if ~strcmpi(displayGeometry, geometry)
   wid = sprintf('%s:%s:ignoringDisplayType', getcomp, fcnName);
   warning(wid, 'Function %s%s%s%s\n%s%s%s\n%s%s\n',  upper(fcnName), ...
      ' expected DisplayType ''', displayGeometry, '''', ...
      'to match Geometry ''', geometry, '''.', ...
      upper(fcnName), ' is ignoring the DisplayType value.');
end
