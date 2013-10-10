function varargout = scaleruler(varargin)
%SCALERULER Add or modify graphic scale on map axes
%
%   SCALERULER toggles the display of a graphic scale.  A graphic scale is
%   a ruler-like graphic object that shows distances on the ground at the
%   correct size for the projection.  If no graphic scale is currently
%   displayed in the current map axes, one will be added.  If any graphic
%   scales are currently displayed, they will be removed.
%
%   SCALERULER ON adds a graphic scale to the current map axes. Multiple
%   graphic scales can be added to the same map axes.
%
%   SCALERULER OFF removes any currently displayed graphic scales.
%
%   SCALERULER('property',value,...) adds a graphic scale and sets the
%   properties to the values specified. A list of graphic scale properties
%   is  displayed by the command SETM(h), where h is the handle to a
%   graphic scale object. The current values for a displayed graphic scale
%   object can be retrieved using GETM.  The properties of a displayed
%   graphic scale object can be modified using SETM.
%
%   Modifying the properties of the graphic scale with SETM results in the
%   replacement of the original object. For this reason, handles to the
%   graphic scale object will change. Use HANDLEM('scaleruler') to get a
%   list of the current handles to graphic scale objects.
%
%   The graphic scale object can be repositioned by dragging the scaleruler
%   with the mouse. The position can also be changed by modifying the
%   'XLoc' and 'YLoc' properties using SETM.
%
%   H = SCALERULER(...) returns the hggroup handle to the graphic scale
%   object.
%
%   See also DISTANCE, SURFDIST, AXESSCALE, PAPERSCALE, DISTORTCALC,
%            MDISTORT.

% Copyright 1996-2009 The MathWorks, Inc.
% $Revision: 1.6.4.10 $  $Date: 2009/03/30 23:39:10 $

if nargin == 2
    s = varargin{1};
    action = varargin{2};

    if ischar(action) && isstruct(s)
        switch lower(action)
            case 'defaults' % issued by the code
                s = defaultscale(s); 
                varargout{1} = s; 
                return
            case 'plot' 	% issued by the code
                h = makescaleruler(s); 
                varargout{1} = h; 
                return
            otherwise 		% property value pairs
                h = scaleruler('on'); 	% create scaleruler
                setm(h,s,action); 	% change the property
                if nargout == 1
                    varargout{1} = h; 
                end
        end
        
    else % property value pairs
        h = scaleruler('on'); 	% create scaleruler
        setm(h,s,action); 	    % change the property
        if nargout == 1; 
            varargout{1} = h; 
        end
        return
    end
elseif nargin == 1 % from the command line, as in 'scaleruler off'

    action = varargin{1};
    switch lower(action)
        case 'off'	% Maybe this should hide, rather than delete

            hndls = [];
            for i=1:20 % don't expect more than 20 distinct scalerulers
                tagstr = ['scaleruler' num2str(i)];
                hexists = findall(gca,'tag',tagstr,'HandleVisibility','on','type','hggroup');
                hndls = [hndls hexists]; %#ok<AGROW>
            end

            for i=1:length(hndls)
                s = get(hndls(i),'Userdata'); % Get the properties structure
                childtag = s.Children; % Handles of all elements of scale ruler
                hchild = findall(gca,'tag',childtag);
                delete([hndls(i); hchild(:)]);
            end
            return
        case 'on'
            % drop down to create scaleruler
        otherwise
            error(['map:' mfilename ':mapdispError'], ...
                'Recognized scaleruler actions are ON and OFF')
    end
elseif nargin == 0

    hndls = [];
    for i=1:20 % don't expect more than 20 distinct scalerulers
        tagstr = ['scaleruler' num2str(i)];
        hexists = findall(gca,'tag',tagstr,'HandleVisibility','on');
        hndls = [hndls hexists]; %#ok<AGROW>
    end

    if isempty(hndls)
        scaleruler('on');
    else
        scaleruler('off');
    end
    return

elseif nargin > 2 && ~mod(nargin,2) % property value pairs, in pairs
    h = scaleruler('on'); 	% create scaleruler
    setm(h,varargin{:}); 	% change the properties
    if nargout == 1; 
        varargout{1} = h; 
    end
    return
    
elseif nargin ~= 0
    error(['map:' mfilename ':mapdispError'], ...
        'Incorrect number of arguments')
end

% check for globe
if strmatch(getm(gca,'mapprojection'),'globe')
    error(['map:' mfilename ':mapdispError'], ...
        'Map scale ruler not appropriate for Globe')
end

% build up structure with properties defining the scale
s = initscale;
s = defaultscale(s);

% display the scale based on the properties
h = makescaleruler(s);

% output argument
if nargout == 1;
    varargout{1} = h;
elseif nargout ~= 0
    error(['map:' mfilename ':mapdispError'], ...
        'Incorrect number of output arguments')
end

%--------------------------------------------------------------------------
function s = initscale

s = struct(...
    'Azimuth',				[],...
    'Children',             [],...
    'Color',				[],...
    'FontAngle',			[],...
    'FontName',				[],...
    'FontSize',				[],...
    'FontUnits',			[],...
    'FontWeight',			[],...
    'Label',				[],...
    'Lat',					[],...
    'Long',					[],...
    'LineWidth',			[],...
    'MajorTick',			[],...
    'MajorTickLabel',		[],...
    'MajorTickLength',		[],...
    'MinorTick',			[],...
    'MinorTickLabel',		[],...
    'MinorTickLength',		[],...
    'Radius',				[],...
    'RulerStyle',			[],...
    'TickDir',				[],...
    'TickMode',				[],...
    'Units',				[],...
    'XLoc',					[],...
    'YLoc',					[], ...
    'ZLoc',					[] ...
    );

%--------------------------------------------------------------------------
function s = defaultscale(s)

if isempty(s.Color);				s.Color = [0 0 0];				end
if isempty(s.FontAngle);			s.FontAngle = 'normal';			end
if isempty(s.FontName);				s.FontName = 'Helvetica';		end
if isempty(s.FontSize);				s.FontSize = 9;				    end
if isempty(s.FontUnits);			s.FontUnits = 'points';			end
if isempty(s.FontWeight);			s.FontWeight = 'normal';		end
if isempty(s.Label);				s.Label = '';					end
if isempty(s.Lat);					s.Lat = [];						end
if isempty(s.LineWidth);			s.LineWidth = 0.5;			    end
if isempty(s.Long);					s.Long = [];					end
if isempty(s.MajorTick);			s.MajorTick = [];				end
if isempty(s.MajorTickLabel);		s.MajorTickLabel = [];			end
if isempty(s.MajorTickLength);		s.MajorTickLength = [];			end
if isempty(s.MinorTick);			s.MinorTick = [];				end
if isempty(s.MinorTickLabel);		s.MinorTickLabel = [];			end
if isempty(s.MinorTickLength);		s.MinorTickLength = [];			end
if isempty(s.Radius);				s.Radius = 'earth';				end
if isempty(s.RulerStyle);			s.RulerStyle = 'ruler';			end
if isempty(s.TickMode);		        s.TickMode = 'auto';			end
if isempty(s.TickDir);		        s.TickDir = 'up';				end
if isempty(s.Units);				s.Units = 'km';					end
if isempty(s.XLoc);					s.XLoc = [];					end
if isempty(s.YLoc);					s.YLoc = [];					end


% Default location at which scale is computed: just off center of x and y limits
if isempty(s.Lat) || isempty(s.Long) || isempty(s.Azimuth)
    xlim = get(gca,'xlim');
    ylim = get(gca,'ylim');
    x = min(xlim) + .4*abs(diff(xlim));
    y = min(ylim) + .4*abs(diff(ylim));
    if isempty(s.Lat) || isempty(s.Long)
        [s.Lat,s.Long] = minvtran(x,y);
    end
    if isempty(s.Azimuth)
        s.Azimuth = 0;
    end
end

% Default tic mark locations
if isempty(s.MajorTick) || isempty(s.MinorTick) || isempty(s.TickMode);
    s.TickMode = 'auto';
end

if strcmp(s.TickMode,'auto')
    nmajticks = 5;
    nminticks = 4;

    angleunits = getm(gca,'angleunits');

    flatlim = getm(gca,'fLatLim');
    flatlim(isinf(flatlim)) = 0;

    flonlim = getm(gca,'fLonLim');

    lonrange = min(abs(diff(flatlim)),abs(diff(flonlim)));

    % convert to degrees and then to distance units
    lonrange = deg2dist(toDegrees(angleunits,lonrange),s.Units);
    
    maxtick = (10^floor(log10(lonrange)));
    MajorTickLim = [ 0 maxtick/2];

    majinc = maxtick/nmajticks;
    MinorTickLim = [ 0 majinc/2];

    % construct major tics and ticlabels in auto mode
    lim = MajorTickLim;
    inc = abs(lim(2))/nmajticks;
    s.MajorTick = lim(1):inc:lim(2);

    % construct minor tics  in auto mode
    lim = MinorTickLim;
    inc = abs(lim(2))/nminticks;
    s.MinorTick = lim(1):inc:lim(2);
    s.MajorTickLabel = [];
end

% Major tick labels
if isempty(s.MajorTickLabel)
    s.MajorTickLabel = num2cell(num2str(s.MajorTick'),2); % need to remove leading blanks
end

% Minor tick labels
if isempty(s.MinorTickLabel)
    tics = s.MinorTick';
    s.MinorTickLabel = num2str(tics(end));
end

% Tick lengths
if isempty(s.MajorTickLength)
    tics = s.MajorTick;
    inc = diff(tics(1:2));
    s.MajorTickLength = inc/5;
end
if isempty(s.MinorTickLength)
    tics = s.MinorTick;
    inc = diff(tics(1:2));
    s.MinorTickLength = inc/2;
end

% default location of scale: 15 percent in from lower left corner of axes
if isempty(s.XLoc) || isempty(s.YLoc)
    xlim = get(gca,'xlim');
    ylim = get(gca,'ylim');
    limx = xlim;
    s.XLoc = limx(1) + 0.15*diff(limx);
    limy = ylim;
    s.YLoc = limy(1) + 0.15*diff(limy);
end

%--------------------------------------------------------------------------
function hout = makescaleruler(s)
% Support multiple styles of the rulers

switch s.RulerStyle
    case 'ruler'   % Construct standard ruler
        hsave = normalruler(s);

    case 'lines'
        hsave = lineruler(s);

    case 'patches'
        hsave = patchruler(s);
end

try
    zdatam(hsave,s.ZLoc)
catch
end

% Construct a unique tag string for the children of the scaleruler baseline
i=0;
while 1
    i=i+1;
    tagstr = ['scaleruler' num2str(i)];
    hexists = findall(gca,'tag',tagstr);
    if isempty(hexists); 
        break; 
    end
end

% save handles
hsave = sort(hsave);
s.Children = tagstr;

set(hsave,'handlevisibility','on','tag',tagstr,'buttondownfcn',@downscale)
hg=hggroup;
set(hsave, 'parent',hg);
set(hg,'userdata',s,'buttondownfcn',@downscale,'hittestarea','on', 'Tag', tagstr);
hout = hg;

%--------------------------------------------------------------------------
function hsave = normalruler(s)
% Determine scaling between geographic and surface units

% Cartesian coordinates of starting point
[xo,yo] = mfwdtran(s.Lat,s.Long);

% Geographical and paper coordinates of a point downrange
sdistdeg = dist2deg(1,s.Units,s.Radius);

[nlat,nlon] = reckon(s.Lat,s.Long,sdistdeg,s.Azimuth);
[xn,yn] = mfwdtran(nlat,nlon); % cartesian coordinates

% Distance between two points in cartesian coordinats
dDdS = sqrt((xo-xn)^2+(yo-yn)^2); % data distance per unit surfacedistance

% Allow for tics and text below the baseline
switch s.TickDir
    case 'down'
        ydir = -1;
        valign = 'Top'; % vertical alignment of text
        labelvalign = 'Bottom';
    otherwise
        ydir = 1;
        valign = 'Bottom';
        labelvalign = 'Top'; % vertical alignment of text
end
hsave = [];

% Minor ticks
xminor = s.XLoc + dDdS*s.MinorTick;
yminor = s.YLoc*ones(size(xminor));

% Construct Major ticks
x = s.XLoc + dDdS*s.MajorTick; %  + dDdS*s.MinorTick(end)
y = s.YLoc*ones(size(x));
y2 = y + ydir*s.MajorTickLength*dDdS;

% Base line
h = plot([xminor x],[yminor y],'k','linewidth',s.LineWidth); hsave = [hsave;h];

% Major ticks
h = plot([x;x],[y;y2],'k','linewidth',s.LineWidth); hsave = [hsave;h];

% Legend label
h = text(mean([x(1) x(end)]), y(1) - ydir*s.MajorTickLength*dDdS, s.Label,...
    'HorizontalAlignment','center','VerticalAlignment',labelvalign,...
    'FontUnits',s.FontUnits,'FontSize',s.FontSize,'Color',s.Color,'FontWeight',s.FontWeight,'FontAngle',s.FontAngle,'FontName',s.FontName);
hsave = [hsave;h];

% Text labels
for i=1:length(s.MajorTickLabel)-1
    h = text(x(i),y2(i),leadblnk(s.MajorTickLabel{i}),...
        'HorizontalAlignment','center','VerticalAlignment',valign,...
        'FontUnits',s.FontUnits,'FontSize',s.FontSize,'FontWeight',s.FontWeight,'FontAngle',s.FontAngle,'FontName',s.FontName);
    hsave = [hsave;h]; %#ok<AGROW>
end

fill=char(32*ones(1,2+length(s.Units)));
h = text(x(end),y2(end),[fill leadblnk(s.MajorTickLabel{end}) ' ' s.Units],...
    'HorizontalAlignment','center','VerticalAlignment',valign,...
    'FontUnits',s.FontUnits,'FontSize',s.FontSize,'FontWeight',s.FontWeight,'FontAngle',s.FontAngle,'FontName',s.FontName);
hsave = [hsave;h];

% Minor ticks
x = s.XLoc - dDdS*s.MinorTick;
y = s.YLoc*ones(size(x));
y2 = y + ydir*s.MinorTickLength*dDdS;
h = plot([x;x],[y;y2],'k','linewidth',s.LineWidth);
hsave = [hsave;h];

% Base line
h = plot(x(1:end),y(1:end),'k','linewidth',s.LineWidth); hsave = [hsave;h];

% Minor tick text label
y2 = y + ydir*s.MajorTickLength*dDdS;
h = text(x(end),y2(end),leadblnk(s.MinorTickLabel),...
    'HorizontalAlignment','center','VerticalAlignment',valign,...
    'FontUnits',s.FontUnits,'FontSize',s.FontSize,'FontWeight',s.FontWeight,'FontAngle',s.FontAngle,'FontName',s.FontName);
hsave = [hsave;h];

% Misc properties
set(hsave,'Color',s.Color);

%--------------------------------------------------------------------------
function hsave = lineruler(s)
% Determine scaling between geographic and surface units

% Cartesian coordinates of starting point
[xo,yo] = mfwdtran(s.Lat,s.Long);

% Geographical and paper coordinates of a point downrange
sdistdeg = dist2deg(1,s.Units,s.Radius);

[nlat,nlon] = reckon(s.Lat,s.Long,sdistdeg,s.Azimuth);
[xn,yn] = mfwdtran(nlat,nlon); % cartesian coordinates

% Distance between two points in cartesian coordinats
dDdS = sqrt((xo-xn)^2+(yo-yn)^2); % data distance per unit surfacedistance

% Allow for tics and text below the baseline
switch s.TickDir
    case 'down'
        ydir = -1;
        valign = 'Top'; % vertical alignment of text
        labelvalign = 'Bottom';
    otherwise
        ydir = 1;
        valign = 'Bottom';
        labelvalign = 'Top'; % vertical alignment of text
end
hsave = [];

% Minor ticks
xminor = s.XLoc - dDdS*s.MinorTick;
yminor = s.YLoc*ones(size(xminor));

% Construct Major ticks
x = s.XLoc + dDdS*s.MajorTick; %  + dDdS*s.MinorTick(end)
y = s.YLoc*ones(size(x));
y2 = y + ydir*s.MajorTickLength*dDdS;

% Base line
h = plot([xminor x],[yminor y],'k','linewidth',s.LineWidth); hsave = [hsave;h];
h = plot(x(1:end),y2(1:end),'k','linewidth',s.LineWidth); hsave = [hsave;h];

% Plot Major ticks
h = plot([x;x],[y;y2],'k','linewidth',s.LineWidth); hsave = [hsave;h];

% Center bars: indices of every other pair
ind =floor(2.5:.5:length(s.MajorTick));
ind = reshape(ind,2,length(ind)/2);
ind = ind(:,1:2:size(ind,2));
if size(ind,2) == 1;
    ind = [ind ind];
end

if ~isempty(ind)
    xbars = x(ind);
    xbars = [xbars ; nan*ones(size(xbars(1,:)))];
    xbars = xbars(:);

    ymid = mean([y;y2]);
    ybars = ymid(ind);
    ybars = [ybars ; nan*ones(size(ybars(1,:)))];
    ybars = ybars(:);

    h = plot(xbars,ybars,'k','linewidth',s.LineWidth); hsave = [hsave;h];
end

% Legend label
h = text(mean([x(1) x(end)]), y(1) - ydir*s.MajorTickLength*dDdS, s.Label,...
    'HorizontalAlignment','center','VerticalAlignment',labelvalign,...
    'FontUnits',s.FontUnits,'FontSize',s.FontSize,'Color',s.Color,'FontWeight',s.FontWeight,'FontAngle',s.FontAngle,'FontName',s.FontName);
hsave = [hsave;h];

% Text labels
for i=1:length(s.MajorTickLabel)-1
    h = text(x(i),y2(i),leadblnk(s.MajorTickLabel{i}),...
        'HorizontalAlignment','center','VerticalAlignment',valign,...
        'FontUnits',s.FontUnits,'FontSize',s.FontSize,'FontWeight',s.FontWeight,'FontAngle',s.FontAngle,'FontName',s.FontName);
    hsave = [hsave;h]; %#ok<AGROW>
end

fill=char(32*ones(1,2+length(s.Units)));
h = text(x(end),y2(end),[fill leadblnk(s.MajorTickLabel{end}) ' ' s.Units],...
    'HorizontalAlignment','center','VerticalAlignment',valign,...
    'FontUnits',s.FontUnits,'FontSize',s.FontSize,'FontWeight',s.FontWeight,'FontAngle',s.FontAngle,'FontName',s.FontName);
hsave = [hsave;h];

% Minor tick center bars: indices of every other pair
x = s.XLoc - dDdS*s.MinorTick;
y = s.YLoc*ones(size(x));
y2 = y + ydir*s.MajorTickLength*dDdS;

ind =floor(1.5:.5:length(s.MinorTick));
ind = reshape(ind,2,length(ind)/2);
ind = ind(:,1:2:size(ind,2));

if size(ind,2) == 1;
    ind = [ind ind];
end

if ~isempty(ind)
    xbars = x(ind);
    xbars = [xbars ; nan*ones(size(xbars(1,:)))];
    xbars = xbars(:);

    ymid = mean([y;y2]);
    ybars = ymid(ind);
    ybars = [ybars ; nan*ones(size(ybars(1,:)))];
    ybars = ybars(:);

    h = plot([x;x],[y;y2],'k','linewidth',s.LineWidth); hsave = [hsave;h];
    h = plot(xbars,ybars,'k','linewidth',s.LineWidth); hsave = [hsave;h];
    h = plot(x(1:end),y(1:end),'k','linewidth',s.LineWidth); hsave = [hsave;h];
    h = plot(x(1:end),y2(1:end),'k','linewidth',s.LineWidth); hsave = [hsave;h];
end

% Minor tick text label
h = text(x(end),y2(end),leadblnk(s.MinorTickLabel),...
    'HorizontalAlignment','center','VerticalAlignment',valign,...
    'FontUnits',s.FontUnits,'FontSize',s.FontSize,'FontWeight',s.FontWeight,'FontAngle',s.FontAngle,'FontName',s.FontName);
hsave = [hsave;h];

% Misc properties
set(hsave,'Color',s.Color);

%--------------------------------------------------------------------------
function hsave = patchruler(s)
% Determine scaling between geographic and surface units

% Cartesian coordinates of starting point
[xo,yo] = mfwdtran(s.Lat,s.Long);

% Geographical and paper coordinates of a point downrange
sdistdeg = dist2deg(1,s.Units,s.Radius);

[nlat,nlon] = reckon(s.Lat,s.Long,sdistdeg,s.Azimuth);
[xn,yn] = mfwdtran(nlat,nlon); % cartesian coordinates

% Distance between two points in cartesian coordinats
dDdS = sqrt((xo-xn)^2+(yo-yn)^2); % data distance per unit surfacedistance

% Allow for tics and text below the baseline
switch s.TickDir
    case 'down'
        ydir = -1;
        valign = 'Top'; % vertical alignment of text
        labelvalign = 'Bottom';
    otherwise
        ydir = 1;
        valign = 'Bottom';
        labelvalign = 'Top'; % vertical alignment of text
end
hsave = [];

% Minor ticks
xminor = s.XLoc - dDdS*s.MinorTick;
yminor = s.YLoc*ones(size(xminor));

% Construct Major ticks
x = s.XLoc + dDdS*s.MajorTick; %  + dDdS*s.MinorTick(end)
y = s.YLoc*ones(size(x));
y2 = y + ydir*s.MajorTickLength*dDdS;

% Base line
h = plot([xminor x],[yminor y],'color',s.Color); hsave = [hsave;h];

h = plot(x(1:end),y2(1:end),'color',s.Color); hsave = [hsave;h];

% Plot Major ticks
h = plot([x;x],[y;y2],'color',s.Color);  hsave = [hsave;h];

% Center bars: indices of every other pair
ind =floor(2.5:.5:length(s.MajorTick));
ind = reshape(ind,2,length(ind)/2);
ind = ind(:,1:2:size(ind,2));

if size(ind,2) == 1;
    ind = [ind ind];
end

xbars = [x(flipud(ind)); x(ind)];
ybars = [y(flipud(ind)); y2(ind)];

hpatch = patch(xbars,ybars,'k'); hsave = [hsave;hpatch];
set(hpatch,'FaceColor',s.Color,'EdgeColor',s.Color);

% Legend label
h = text(mean([x(1) x(end)]), y(1) - ydir*s.MajorTickLength*dDdS, s.Label,...
    'HorizontalAlignment','center','VerticalAlignment',labelvalign,...
    'FontUnits',s.FontUnits,'FontSize',s.FontSize,'Color',s.Color, ...
    'FontWeight',s.FontWeight,'FontAngle',s.FontAngle,'FontName',s.FontName);
hsave = [hsave;h];

% Text labels
for i=1:length(s.MajorTickLabel)-1
    h = text(x(i),y2(i),leadblnk(s.MajorTickLabel{i}),...
        'HorizontalAlignment','center','VerticalAlignment',valign,...
        'FontUnits',s.FontUnits,'FontSize',s.FontSize,'Color',s.Color, ...
        'FontWeight',s.FontWeight,'FontAngle',s.FontAngle,'FontName',s.FontName);
    hsave = [hsave;h]; %#ok<AGROW>
end

fill=char(32*ones(1,2+length(s.Units)));
h = text(x(end),y2(end),[fill leadblnk(s.MajorTickLabel{end}) ' ' s.Units],...
    'HorizontalAlignment','center','VerticalAlignment',valign,...
    'FontUnits',s.FontUnits,'FontSize',s.FontSize,'Color',s.Color, ...
    'FontWeight',s.FontWeight,'FontAngle',s.FontAngle,'FontName',s.FontName);
hsave = [hsave;h];


% Minor tick center bars: indices of every other pair
x = s.XLoc - dDdS*s.MinorTick;
y = s.YLoc*ones(size(x));
y2 = y + ydir*s.MajorTickLength*dDdS;

ind =floor(1.5:.5:length(s.MinorTick));
ind = reshape(ind,2,length(ind)/2);
ind = ind(:,1:2:size(ind,2));

if size(ind,2) == 1;
    ind = [ind ind];
end

xbars = [x(flipud(ind));  x(ind)];
ybars = [y(flipud(ind)); y2(ind)];

hpatch = patch(xbars,ybars,'k'); hsave = [hsave;hpatch];
set(hpatch,'FaceColor',s.Color);

h = plot([x;x],[y;y2],'Color',s.Color); hsave = [hsave;h];
h = plot(x(1:end),y(1:end),'Color',s.Color); hsave = [hsave;h];
h = plot(x(1:end),y2(1:end),'Color',s.Color); hsave = [hsave;h];

% Minor tick text label
h = text(x(end),y2(end),leadblnk(s.MinorTickLabel),...
    'HorizontalAlignment','center','VerticalAlignment',valign,...
    'FontUnits',s.FontUnits,'FontSize',s.FontSize,'Color',s.Color,...
    'FontWeight',s.FontWeight,'FontAngle',s.FontAngle,'FontName',s.FontName);
hsave = [hsave;h];

%--------------------------------------------------------------------------
function downscale(action,delta,direction) %#ok<INUSD,INUSD>
%  Mouse clicks on the scale baseline

%  Ensure that the plot is in a 2D view
if any(get(gca,'view') ~= [0 90])
    btn = questdlg({'Must be in 2D view for operation.',...
        'Change to 2D view?'},...
        'Incorrect View','Change','Cancel','Change');

    switch btn
        case 'Change',      view(2);
        case 'Cancel',      return
    end
end

set(gcf,'WindowButtonMotionFcn',@movescale,...
    'WindowButtonUpFcn',@upscale,...
    'pointer','fleur')

hobject = getScaleRulerObject;

pt = get(gca,'CurrentPoint');
s  = get(hobject,'UserData');
s.oldUD = s;
s.pt    = pt;

set(hobject, 'UserData', s);
set(hobject, 'selected', 'on');
set(hobject, 'Visible',   'on');

%--------------------------------------------------------------------------
function movescale(action,delta,direction) %#ok<INUSD,INUSD>
%  Move the origin marker by the mouse

hobject = getScaleRulerObject;

% Update the location of the scaleruler.
% The location is set with XLoc and YLoc.
t = get(hobject,'UserData');
lastpt = t.pt;
pt = get(gca,'CurrentPoint');
t.pt = pt;
delx = pt(1,1) - lastpt(1,1);
dely = pt(1,2) - lastpt(1,2);
t.XLoc = t.XLoc + delx;
t.YLoc = t.YLoc + dely;
set(hobject,'UserData',t);

children = get(hobject,'Children');

% Line
lines = findobj(children,'Type','line');
resetXYData(lines, delx, dely);

% Patch
patches = findobj(children,'Type','patch');
resetXYData(patches, delx, dely);

% Text
textObjects = findobj(children,'Type','text');
for j=1:numel(textObjects);
    position = get(textObjects(j),'Position') +[delx, dely, 0.0];
    set(textObjects(j),'Position',position);
end

%--------------------------------------------------------------------------
function resetXYData(obj, deltaX, deltaY)
% Reset the X and Y Data by the values deltaX and deltaY

for j=1:numel(obj);
    x  = get(obj(j),'Xdata');
    y  = get(obj(j),'Ydata');
    set(obj(j), 'Xdata', x+deltaX, 'Ydata', y+deltaY);
end

%--------------------------------------------------------------------------
function upscale(action,delta,direction) %#ok<INUSD,INUSD>
%  Mouse release of the origin marker

hobject   = getScaleRulerObject;
set(gcf,'WindowButtonMotionFcn','',...
    'WindowButtonUpFcn','',...
    'pointer','arrow')
set(hobject,'selected','off')

%--------------------------------------------------------------------------
function handle = getScaleRulerObject
tag = get(gco,'Tag');
handle = findobj(gca, 'Tag', tag, 'Type', 'hggroup');

%--------------------------------------------------------------------------
function dist = deg2dist(deg,units)
% Convert from spherical distance in degrees
%
%   DIST = DEG2DIST(DEG, UNITS) converts distances from degrees, as
%   measured along a great circle on a sphere with a radius of 6371 km
%   (the mean radius of the Earth) to the UNITS of length or angle
%   specified by the string  in UNITS.  If UNITS is 'degrees' or
%   'radians', then DIST itself will be a spherical distance.

angleUnits = {'degrees','radians'};
k = find(strncmpi(deblank(units), angleUnits, numel(deblank(units))));
if numel(k) == 1
    % In case units is 'degrees' or 'radians'
    dist = fromDegrees(angleUnits{k}, deg);
else
    % Assume that units specifies a length unit; convert using
    % kilometers as an intermediate unit.
    dist = unitsratio(units,'km') * deg2km(deg);
end
