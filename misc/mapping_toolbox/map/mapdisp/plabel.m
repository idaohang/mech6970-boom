function hndl = plabel(varargin)
%PLABEL Toggle and control display of parallel labels
%
%   PLABEL toggles the display of the parallel labels on the map axes.
%   These labels are drawn using the properties specified in the map axes.
%
%   PLABEL ON turns the parallel labels on. PLABEL OFF turns them off.
%
%   PLABEL RESET will redraw the parallel labels with the currently
%   specified properties.  This differs from the ON and OFF which simply
%   sets the visible property of the current labels.
%
%   PLABEL(meridian) places the parallel labels at the specified meridian.
%   The input meridian is used to set the PLabelMeridian property in the
%   map axes.
%
%   PLABEL('MapAxesPropertyName',PropertyValue,...) uses the specified Map
%   Axes properties to draw the parallel labels.
%
%   H = PLABEL(...) returns the handles of the labels drawn.
%
%   See also MLABEL, AXESM, SETM, SET.

% Copyright 1996-2009 The MathWorks, Inc.
% $Revision: 1.10.4.7 $  $Date: 2009/03/09 19:16:34 $
% Written by:  E. Byrns, E. Brown

mstruct = gcm;

h = handlem('PLabel');
if nargout ~= 0
    hndl = h;
end

if nargin == 0
    if ~isempty(h)
	    if strcmp('off',get(h,'Visible'))
	          showm('PLabel');
			  mstruct.parallellabel = 'on';
              set(gca,'UserData',mstruct)
			  return
	    else
	          hidem('PLabel');
			  mstruct.parallellabel = 'off';
              set(gca,'UserData',mstruct)
			  return
	    end
    end

elseif nargin == 1 && strcmpi(varargin{1},'on')
    if ~isempty(h)                      %  Show existing parallel labels.
 	      showm('PLabel');               %  Else, draw new one
		  mstruct.parallellabel = 'on';
		  set(gca,'UserData',mstruct)
		  return
    end

elseif nargin == 1 && strcmpi(varargin{1},'off')
 	hidem('PLabel');
    mstruct.parallellabel = 'off';
    set(gca,'UserData',mstruct)
    return

elseif nargin == 1 && ~strcmpi(varargin{1},'reset')
    % AXESM recursively calls PLABEL to display the labels
    axesm(mstruct,'ParallelLabel','reset','PLabelMeridian',varargin{1});
    return

elseif rem(nargin,2) == 0
    % AXESM recursively calls PLABEL to display the labels
    axesm(mstruct,'ParallelLabel','reset',varargin{:});
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
    case 'compass',   format = 'ns';
	case 'signed',    format = 'pm';
	otherwise,        format = 'none';
end

%  Get the parallel label properties

pposit  = mstruct.plabellocation;
pplace  = mstruct.plabelmeridian;
pround  = mstruct.plabelround;

%  Get the necessary current map data

maplat  = mstruct.maplatlimit;
origin  = mstruct.origin;
units   = mstruct.angleunits;
frmlon  = mstruct.flonlimit;
gridalt = mstruct.galtitude;

%  Set grid to above top of z axis if altitude is set to inf.

if isinf(gridalt);   gridalt = max(get(gca,'Zlim'))+1;   end

%  Convert the input data into degrees.
%  DMS presents problems with arithmetic below

[maplat,origin,frmlon,pposit,pplace] = toDegrees(units,maplat,origin,frmlon,pposit,pplace);
epsilon = 500*epsm('degrees');

%  Skip labeling if inf or NaN entered

if any(isinf(pposit)) || any(isnan(pposit))
    if nargout == 1
        hndl = [];
    end
    return
end

%  Latitude locations for the whole world

latlim = [-90 90];

%  Compute the latitudes at which to place labels

if length(pposit) == 1
    latline = [fliplr(-pposit:-pposit:min(latlim)), 0:pposit:max(latlim) ];
else
	latline = pposit;            %  Vector of points supplied
end

latline = latline(latline >= min(maplat) & latline <= max(maplat));

%  Compute the latitude placement points

lonline = pplace(ones(size(latline)));

%  Set appropriate horizontal justification

if pplace == min(frmlon) + origin(2)
     justify = 'right';
     if  min(frmlon) == -180
	     lonline = lonline + epsilon;   %  Slightly inside the frame
     else
         lonline = lonline - epsilon;   %  Slightly outside the frame
     end

elseif pplace == max(frmlon) + origin(2)
     justify = 'left';
     if max(frmlon) == 180
         lonline = lonline - epsilon;   %  Slightly inside the frame
	 else
         lonline = lonline + epsilon;   %  Slightly outside the frame
     end
else
     justify = 'center';
end

%  Compute the label string matrix

parallels = latline;
if strncmpi(labelunits, 'dms', numel(labelunits))
    % Replace 'dms' with 'degrees2dms'
    % and     'dm'  with 'degrees2dm'
    labelunits = ['degrees2' labelunits];
else
    % labelunits should be 'degrees' or 'radians'
    parallels = fromDegrees(labelunits,parallels);
end
labelstr = angl2str(parallels,format,labelunits,pround);

%  Transform the location data back into the map units

[latline,lonline] = fromDegrees(units,latline,lonline);

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
               'HorizontalAlignment',justify,...
               'Interpreter','tex',...
		       'VerticalAlignment','middle',...
		       'Tag','PLabel',...
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

mstruct.parallellabel = 'on';  set(gca,'UserData',mstruct)

%  Set handle return argument if necessary

if nargout == 1
    hndl = hndl0;
end
