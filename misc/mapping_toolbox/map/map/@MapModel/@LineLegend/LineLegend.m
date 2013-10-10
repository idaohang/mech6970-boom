function h = LineLegend(varargin)
%LINELEGEND Define symbolization for line data.
%
%   L = LINELEGEND constructs a linelegend with default values.
%
%   L =
%   LINELEGEND('PropertyName1',PropertyValue1,'PropertyName2',PropertyValue2,...)
%   Constructs a linelegend and sets multiple property values.  Valid
%   PropertyNames are Color, LineStyle, and LineWidth. PropertyValue may be a
%   function handle.  The function must have one input, typically an attribute
%   value, and the function must return a valid property value for the
%   associated PropertyName.

% Copyright 1996-2007 The MathWorks, Inc.
% $Revision: 1.1.6.4 $  $Date: 2007/11/09 20:25:59 $

h = MapModel.LineLegend;

% Default Line Properties
h.Color     = {'Default','',h.getColor};
h.LineStyle = {'Default','','-'};
h.LineWidth = {'Default','',0.5};
h.Visible   = {'Default','','on'};

if nargin > 0
  % User defined properties
  fldnames = fieldnames(varargin{1});
  for i=1:length(fldnames)
    if isprop(h,fldnames{i})
      h.(fldnames{i}) = varargin{1}.(fldnames{i});
    else
      error(['map:' mfilename ':mapError'], ...
          '%s is not a property that can be set for a Line shape.', ...
          fldnames{i})
    end
  end
end