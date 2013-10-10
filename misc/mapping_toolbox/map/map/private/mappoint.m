function varargout = mappoint(xdata, ydata, varargin)
%MAPPOINT Display points without projection
%
%   Construct a line object for display in map (x-y) coordinates. Set the
%   LineStyle to 'none' so that only the coordinate points are displayed.
%   Set the Marker and MarkerEdgeColor to default values that can be reset
%   with VARARGIN.  Return empty if XDATA and YDATA are empty.
%
%   Example
%   -------
%   load coast
%   figure
%   mappoint(long, lat)

% Copyright 2006 The MathWorks, Inc.
% $Revision: 1.1.6.1 $  $Date: 2006/05/24 03:35:42 $

h = mapline(xdata, ydata, ...
   'Marker', '+', 'MarkerEdgeColor', 'red', varargin{:}, 'LineStyle', 'none');

% Suppress output if called with no return value and no semicolon.
if nargout > 0
    varargout{1} = h;
end
