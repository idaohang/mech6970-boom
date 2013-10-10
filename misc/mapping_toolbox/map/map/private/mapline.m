function varargout = mapline(xdata, ydata, varargin)
%MAPLINE Display line without projection
%
%   Construct a line object for display in map (x-y) coordinates. Return
%   empty if XDATA and YDATA are empty.
%
%   Example
%   -------
%   load coast
%   figure
%   mapline(long, lat, 'Color', 'black')

% Copyright 2006 The MathWorks, Inc.
% $Revision: 1.1.6.1 $  $Date: 2006/05/24 03:35:41 $

% Verify NaN locations are equal.
if ~isequal(isnan(xdata), isnan(ydata))
   eid = sprintf('%s:%s:inconsistentXY', getcomp, mfilename);
   error(eid,'XDATA and YDATA mismatch in size or NaN locations.')
end

if ~isempty(xdata) || ~isempty(ydata)
   % Create the line object and set the default color to blue.
   h = line(xdata(:), ydata(:), 'Color', [0 0 1], varargin{:});
else
   % Either xdata or ydata are empty.
   h = reshape([],[0 1]);
end

% Suppress output if called with no return value and no semicolon.
if nargout > 0
   varargout{1} = h;
end
