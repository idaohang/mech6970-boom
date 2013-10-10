function [h,msg] = linem(varargin)
%LINEM Project line object on map axes
%
%  LINEM(LAT,LON) projects the line objects onto the current map axes.
%  The input latitude and longitude data must be in the same units as
%  specified in the current map axes.  Unlike PLOTM and PLOT3M, LINEM
%  will always add lines to the current axes, regardless of the current
%  hold state.
%
%  LINEM(LAT,LON,Z) projects a 3-D line object onto the current map
%  axes.  The units of Z are arbitrary, except when using the 'globe'
%  projection.  In the case of 'globe', Z should have the same units as
%  the radius of the earth or semimajor axis specified in the 'geoid'
%  (reference ellipsoid) property of the map axes.  This implies that
%  for a reference ellipsoid vector of [1 0] (a unit sphere), the units
%  of Z are earth radii.
%
%  LINEM(LAT,LON,'LineSpec') and LINEM(LAT,LON,Z,'LineSpec') uses any
%  valid LineSpec string to display the line object.
%
%  LINEM(LAT,LON,'PropertyName',PropertyValue,...) and
%  LINEM(LAT,LON,Z,'PropertyName',PropertyValue,...) uses the line
%  object properties specified to display the line objects.  Except for
%  xdata, ydata and zdata, all line properties, and styles available
%  through LINE are supported by LINEM.
%
%  h = LINEM(...) returns the handles to the line objects displayed.
%
%  See also PLOTM, PLOT3M, LINE.

% Copyright 1996-2009 The MathWorks, Inc.
% $Revision: 1.15.4.7 $  $Date: 2009/05/14 17:06:20 $

% Obsolete syntax
% ---------------
% [h,msg] = LINEM(...) returns a string indicating any error encountered.
if nargout > 1
    warnObsoleteMSGSyntax(mfilename)
    msg = '';
end

if nargin == 0
    errorOnMissingUI(mfilename)
end

error(nargchk(2,inf,nargin, 'struct'))

lat = varargin{1};
lon = varargin{2};
if nargin == 2 || ischar(varargin{3})
    z = zeros(size(lat));
    varargin(1:2) = [];
else
    z = varargin{3};
    varargin(1:3) = [];
end

%  Test for scalar z data
if numel(z) == 1;
    z = z + zeros(size(lat));
end

%  Argument size tests
if any([ndims(lat) ndims(lon) ndims(z)] > 2)
    error(['map:' mfilename ':nDimsExceedsTwo'],...
        'Line data must not contain pages.')
elseif ~isequal(size(lat),size(lon),size(z))
    error(['map:' mfilename ':inconsistentDims'],...
        'Inconsistent dimensions on input data.');
end

%  Ensure a column vector if a row vector is given.
if size(lat,1) == 1;
    lat = lat(:);
    lon = lon(:);
    z = z(:);
end

%  Parse the line styles
if rem(length(varargin),2)
    [lstyle,lcolor,lmark,err] = colstyle(varargin{1});
    error(err)  %#ok<ERTAG>
 
    if isempty(lstyle);
        lstyle = 'none';
    end

    if isempty(lmark);
        lmark = 'none';
    end

    %  Build up a new property string vector
    varargin(1) = [];
    linespec = {};
    if ~isempty(lcolor)
        linespec{length(linespec)+1} = 'Color';
        linespec{length(linespec)+1} = lcolor;
    end

    if ~strcmp(lstyle,'none') || ~strcmp(lmark,'none')
        linespec{length(linespec)+1} = 'LineStyle';
        linespec{length(linespec)+1} = lstyle;
        linespec{length(linespec)+1} = 'Marker';
        linespec{length(linespec)+1} = lmark;
    end

    %  Append linespec to front of varargin.  Allows users to
    %  override linespec properties
    varargin = [linespec varargin];
end

%  Validate map axes, project lines, and display
ax = getParentAxesFromArgList(varargin);
if ~isempty(ax)
    mstruct = gcm(ax);
    [x,y,z,savepts] = mfwdtran(mstruct,lat,lon,z,'linem');
    h0 = line(x,y,z,'ButtonDownFcn','uimaptbx','Parent',ax);
else
    mstruct = gcm;
    [x,y,z,savepts] = mfwdtran(mstruct,lat,lon,z,'linem');
    h0 = line(x,y,z,'ButtonDownFcn','uimaptbx');
end

%  Set the userdata property for each line
for i = 1:size(x,2)
    if isempty(savepts.clipped)
        userdata.clipped = [];
    else
        userdata.clipped ...
            = savepts.clipped((savepts.clipped(:,2) == i),:);
    end

    if isempty(savepts.trimmed)
        userdata.trimmed = [];
    else
        userdata.trimmed ...
            = savepts.trimmed((savepts.trimmed(:,2) == i),[1 3 4]);
    end

    set(h0(i),'UserData',userdata)
end

%  Set line properties if necessary
if ~isempty(varargin)
    set(h0,varargin{:});
end

% Assign output arguments if specified
if nargout > 0
    h = h0;
end
