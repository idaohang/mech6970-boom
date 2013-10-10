function [h,msg] = textm(varargin)
%TEXTM Project text annotation on map axes
%
%   TEXTM(LAT, LON, STRING) projects the text in STRING onto the current
%   map axes at the locations specified by the LAT and LON.  The units of
%   LAT and LON must match the 'angleunits' property of the map axes.  If
%   LAT and LON contain multiple elements, TEXTM places a text object at
%   each location.  In this case STRING may be a cell array of strings with
%   the same number of elements as LAT and LON.  (For backward
%   compatibility, STRING may also be a 2-D character array such that
%   size(STRING,1) matches numel(LAT)).
%
%   TEXTM(LAT, LON, Z, STRING) draws the text at the altitude(s) specified
%   in Z, which must be the same size as LAT and LON.
%
%   TEXTM(..., PROP1, VAL1, PROP2, VAL2,...) applies additional text
%   object properties.
%
%   H = TEXTM(...) returns the handles to the text objects drawn.
%
%   Note
%   ----
%   You may be working with scalar LAT and LON data or vector LAT and LON
%   data. If you are in scalar mode and you enter a cell array of strings,
%   you will get a text object with a multiline string. Also note that
%   vertical slash characters, rather than producing multiline strings,
%   will yield a single line string containing vertical slashes. On the
%   other hand, if LAT and LON are nonscalar, then the size of the cell
%   array input must match their size exactly.
%
%   See also GTEXTM, TEXT.

% Copyright 1996-2009 The MathWorks, Inc.
% $Revision: 1.16.4.8 $  $Date: 2009/05/14 17:06:31 $

% Obsolete syntax
% ---------------
%  [H, MSG] = TEXTM(...) returns an optional second output which contains
%  a string indicating any errors encountered.
if nargout > 1
    warnObsoleteMSGSyntax(mfilename)
    msg = '';
end

if nargin == 0
    errorOnMissingUI(mfilename)
end

% Parse input arguments
error(nargchk(3,Inf,nargin,'struct'))

lat = varargin{1};
lon = varargin{2};
if ischar(varargin{3}) || iscell(varargin{3})
    z = zeros(size(lat));
    string = varargin{3};
    propertyValuePairs = varargin(4:end);
elseif nargin >= 4
    z = varargin{3};
    string = varargin{4};
    propertyValuePairs = varargin(5:end);
else
    error(['map:' mfilename ':missingStringArg'], ...
        'Input argument STRING must be provided.')
end

% Check for size consistency
if ~isequal(size(lat),size(lon),size(z))
    error(['map:' mfilename ':coordSizeMismatch'], ...
        'Sizes of LAT, LON, and Z arguments much match.')
end

% From reference page for TEXT: "When specifying the string for a single
% text object, cell arrays of strings and padded string matrices result
% in a text object with a multiline string." This means, for example,
% that it's OK for STRING to be a multi-element cell array if LAT and
% LON are scalar, but their sizes need to match if LAT and LON are
% nonscalar.
if ~isscalar(lat)
    if (iscell(string) && numel(string) ~= numel(lat)) ...
            || (ischar(string) && size(string,1) ~= numel(lat))
        error(['map:' mfilename ':stringSizeMismatch'], ...
            'Size of argument STRING is inconsistent with sizes of coordinate arrays.')
    end
end

% Convert to column vectors
lat = lat(:); 
lon = lon(:);
z = z(:);
if iscell(string)
    string = string(:);
end

%  Validate map axes, project lines, and display
ax = getParentAxesFromArgList(varargin);
if ~isempty(ax)
    mstruct = gcm(ax);
else
    mstruct = gcm;
end

[x,y,z,savepts] = mfwdtran(mstruct, lat, lon, z, 'text');
h0 = text(x, y, z, string,...
    'ButtonDownFcn', 'uimaptbx', 'Clipping', 'on', propertyValuePairs{:});

% Set the userdata property for each line
for i = 1:length(x)
     if isempty(savepts.clipped)
	     userdata.clipped = [];
	 else
	     cindx = (savepts.clipped(:,2) == i);
         userdata.clipped = savepts.clipped(cindx,:);
     end

     if isempty(savepts.trimmed)
	     userdata.trimmed = [];
	 else
	     tindx = (savepts.trimmed(:,2) == i);
         userdata.trimmed = savepts.trimmed(tindx,[1 3 4]);
     end

     set(h0(i),'UserData',userdata)
end

% Assign output arguments if specified
if nargout > 0
    h = h0;
end
