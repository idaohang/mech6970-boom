function hndl = mlabel(varargin)
%MLABEL Toggle and control display of meridian labels
%
%   MLABEL toggles the display of the meridian labels on the map axes.
%   These labels are drawn using the properties specified in the map axes.
%
%   MLABEL ON turns the meridian labels on. MLABEL OFF turns them off.
%
%   MLABEL RESET will redraw the meridian labels with the currently
%   specified properties.  This differs from the ON and OFF option which
%   simply sets the visible property of the current labels.
%
%   MLABEL(parallel) places the meridian labels at the specified parallel.
%   The input parallel is used to set the MLabelParallel property in the
%   map axes.
%
%   MLABEL('MapAxesPropertyName',PropertyValue,...) uses the specified Map
%   Axes properties to draw the meridian labels.
%
%   H = MLABEL(...) returns the handles of the labels drawn.
%
%   See also AXESM, MLABELZERO22PI, PLABEL, SET, SETM.

% Copyright 1996-2009 The MathWorks, Inc.
% $Revision: 1.10.4.7 $  $Date: 2009/03/09 19:16:30 $
% Written by:  E. Byrns, E. Brown

mstruct = gcm;

h = handlem('MLabel');
if nargout ~= 0
    hndl = h;
end

if nargin == 0
    if ~isempty(h)
        if strcmp('off',get(h,'Visible'))
            showm('MLabel');
            mstruct.meridianlabel = 'on';
            set(gca,'UserData',mstruct)
            return
        else
            hidem('MLabel');
            mstruct.meridianlabel = 'off';
            set(gca,'UserData',mstruct)
            return
        end
    end

elseif nargin == 1 && strcmpi(varargin{1},'on')
    if ~isempty(h)                      %  Show existing meridian labels.
 	      showm('MLabel');               %  Else, draw new one
		  mstruct.meridianlabel = 'on';
		  set(gca,'UserData',mstruct)
		  return
    end

elseif nargin == 1 && strcmpi(varargin{1},'off')
 	hidem('MLabel');
    mstruct.meridianlabel = 'off';
    set(gca,'UserData',mstruct)
    return

elseif nargin == 1 && ~strcmpi(varargin{1},'reset')
    % AXESM recursively calls MLABEL to display the labels
    axesm(mstruct,'MeridianLabel','reset','MLabelParallel',varargin{1});
	return        

elseif rem(nargin,2) == 0
    % AXESM recursively calls MLABEL to display the labels
    axesm(mstruct,'MeridianLabel','reset',varargin{:});
    return        

elseif (nargin == 1 && ~strcmpi(varargin{1},'reset') ) || ...
       (nargin > 1 && rem(nargin,2) ~= 0)
    error(['map:' mfilename ':invalidArgCount'], ...
        'Incorrect number of arguments')
end


%  Default operation is to label the map.  Action string = 'reset'

%  Clear existing labels
if ~isempty(h)
    delete(h)
end       

%  Get the font definition properties

fontangle  = mstruct.fontangle;
fontname   = mstruct.fontname;
fontsize   = mstruct.fontsize;
fontunits  = mstruct.fontunits;
fontweight = mstruct.fontweight;
fontcolor  = mstruct.fontcolor;
labelunits = mstruct.labelunits;

%  Convert the format into a string recognized by angl2str

switch mstruct.labelformat
    case 'compass',   format = 'ew';
	case 'signed',    format = 'pm';
	otherwise,        format = 'none';
end

%  Get the meridian label properties

mposit  = mstruct.mlabellocation;
mplace  = mstruct.mlabelparallel;
mround  = mstruct.mlabelround;

%  Get the necessary current map data

maplon  = mstruct.maplonlimit;
units   = mstruct.angleunits;
frmlat  = mstruct.flatlimit;
gridalt = mstruct.galtitude;

%  Set grid to above top of z axis if altitude is set to inf.

if isinf(gridalt);   gridalt = max(get(gca,'Zlim'))+1;   end

%  Convert the input data into degrees.
%  DMS presents problems with arithmetic below

[maplon,frmlat,mposit,mplace] = toDegrees(units,maplon,frmlat,mposit,mplace);
epsilon = 500*epsm('degrees');

%  Skip labeling if inf or NaN entered

if any(isinf(mposit)) || any(isnan(mposit))
    if nargout == 1
        hndl = [];
    end
    return
end

%  Longitude locations for the whole world and then some
%  Will be truncated later.  Use more than the whole world
%  to ensure capturing the data range of the current map.

lonlim = [-360 360];

%  Compute the longitudes at which to place labels

if length(mposit) == 1
	lonline = [fliplr(-mposit:-mposit:min(lonlim)), 0:mposit:max(lonlim) ];
else
	lonline = mposit;            %  Vector of points supplied
end

lonline = lonline(lonline >= min(maplon)  &  lonline <= max(maplon));

%  Compute the latitude placement points and set vertical justification

latline = mplace(ones(size(lonline)));
if mplace == min(frmlat)
     justify = 'top';
     if min(frmlat) == -90
         latline = latline + epsilon;   %  Slightly above south pole
     else
         latline = latline - epsilon;   %  Slightly below bottom of map
     end
elseif mplace == max(frmlat)
     justify = 'bottom';
     if max(frmlat) == 90
         latline = latline - epsilon;   %  Slightly below north pole
     else
         latline = latline + epsilon;   %  Slightly above top of map
     end
else
     justify = 'middle';
end

%  Compute the label string matrix

meridians = npi2pi(lonline,'degrees','exact');
if strncmpi(labelunits, 'dms', numel(labelunits))
    % Replace 'dms' with 'degrees2dms'
    % and     'dm'  with 'degrees2dm'
    labelunits = ['degrees2' labelunits];
else
    % labelunits should be 'degrees' or 'radians'
    meridians = fromDegrees(labelunits,meridians);
end
labelstr = angl2str(meridians,format,labelunits,mround);

%  Transform the location data back into the map units

latline = fromDegrees(units,latline);   %  Avoid a reset of -180 to
lonline = npi2pi(lonline,'degrees','inward');  % +180 which can yield a
lonline = fromDegrees(units,lonline);   %  double label at +180 spot

%  Display the latitude labels on the map

gridalt = gridalt(ones(size(latline)));

if ~isempty(latline) && ~isempty(lonline)

	hndl0 = textm(latline,lonline,gridalt,labelstr,...
               'Color',fontcolor,...
               'FontAngle',fontangle,...
               'FontName',fontname,...
               'FontSize',fontsize,...
               'FontUnits',fontunits,...
               'FontWeight',fontweight,...
               'HorizontalAlignment','center',...
               'Interpreter','tex',...
		       'VerticalAlignment',justify,...
		       'Tag','MLabel',...
			   'Clipping','off');

	%  Align text to graticule
	
	switch mstruct.labelrotation
	case 'on'
		rotatetext(hndl0)
	end

else

	hndl0 = [];

end

%  Set the display flag to on

mstruct.meridianlabel = 'on';  set(gca,'UserData',mstruct)

%  Set handle return argument if necessary

if nargout == 1
    hndl = hndl0;
end
