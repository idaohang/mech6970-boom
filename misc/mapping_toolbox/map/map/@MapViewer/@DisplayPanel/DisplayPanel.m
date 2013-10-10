function this = DisplayPanel(viewer)
%DisplayPanel

% Copyright 2003-2010 The MathWorks, Inc.
% $Revision: 1.1.6.8 $  $Date: 2010/03/22 03:51:52 $

viewerPosition = viewer.getPositionInPixels;

panelHeight = 55;
borderWidth = 2;
pixWidth = viewerPosition(3)-(2*borderWidth);
pixHeight = panelHeight - (2*borderWidth);

x_pos = borderWidth;
y_pos = borderWidth;

this = MapViewer.DisplayPanel;

this.hPanel = uipanel(...
    'Parent',viewer.Figure,...
    'Units','pixel',...
    'BorderType','beveledin',...
    'BorderWidth',borderWidth,...
    'Position',[1, 1, viewerPosition(3), panelHeight],...
    'Visible','on',...
    'HandleVisibility','off', ...
    'DeleteFcn', @deleteUIPanel);

    function deleteUIPanel(hSrc,event) %#ok<INUSD>
        % The UIPanel with handle hPanel will be deleted automatically
        % when its parent figure goes away.  Make sure that the
        % associated DisplayPanel object gets deleted at the same time.
        if ishandle(this)  % Applying ishandle to a non-HG (UDD) object
            delete(this)
        end
    end

this.LayoutPanel = uipanel('Parent',this.hPanel,...
                           'Units','pixels',...
                           'BorderType','none',...
                           'Position',[x_pos, y_pos, pixWidth-2, pixHeight-2],...
                           'Visible','on',...
                           'HandleVisibility','off');

set(viewer.Figure,'Color',get(this.hPanel,'BackgroundColor'));

this.MapUnits = {'None','none';
                 'Kilometers','km';
                 'Meters','m';
%                 'Centimeters','cm';
%                 'Milimeters','mm';
%                 'Microns','u';
                 'Nautical Miles','nm';
                 'International Feet','ft';
%                 'Inches','in';
%                 'Yards','yd';
%                 'International Miles','mi';
                 'US Survey Feet','sf';
%                 'US Survey Miles','sm';
                };

createXYDisplay(this,viewer);
createScaleDisplay(this,viewer);
createMapUnitsDisplay(this,viewer);
createActiveLayerDisplay(this,viewer);
end

%-----------------------------------------------------------------------

function createXYDisplay(this,viewer)

c = get(viewer.Figure,'Color');

pX = [this.PanelXMargin+10,...
      this.PanelYMargin+this.TextBoxHeight-3,...
      this.TextBoxWidth-50,...
      this.TextBoxHeight];

pY = [this.PanelXMargin+10,...
      this.PanelYMargin-3,...
      this.TextBoxWidth-50,...
      this.TextBoxHeight];

uicontrol('Parent',this.LayoutPanel,'Style','text','String','X:',...
          'Units','Pixels','Enable','inactive',...
          'Position',pX+[30 0 0 0],'HorizontalAlignment','left','BackgroundColor',c);

uicontrol('Parent',this.LayoutPanel,'Style','text','String','Y:',...
          'Units','Pixels','Enable','inactive',...
          'Position',pY+[30 0 0 0],'HorizontalAlignment','left','BackgroundColor',c);

p = [this.PanelXMargin+this.TextBoxWidth-50,...
     this.PanelYMargin+this.TextBoxHeight+2,...
     this.TextBoxWidth+40,...
     this.TextBoxHeight];
     
% $$$ p = [this.PanelXMargin + this.TextBoxWidth + spacing,...
% $$$      this.PanelYMargin + this.TextBoxHeight + 2,...
% $$$      this.TextBoxWidth,...
% $$$      this.TextBoxHeight];
viewer.XDisplay = uicontrol(this.LayoutPanel,'Style','edit','Units','Pixels',...
    'String',' ','Enable','inactive','Position',p+[30 0 0 0]);

p = [this.PanelXMargin+this.TextBoxWidth-50,...
     this.PanelYMargin,...
     this.TextBoxWidth+40,...
     this.TextBoxHeight];
     
% $$$ p = [this.PanelXMargin + this.TextBoxWidth + spacing,...
% $$$      this.PanelYMargin,...
% $$$      this.TextBoxWidth,...
% $$$      this.TextBoxHeight];
viewer.YDisplay = uicontrol(this.LayoutPanel,'Style','edit','Units','Pixels',...
    'String',' ','Enable','inactive','Position',p+[30 0 0 0]);
end

%-----------------------------------------------------------------------

function createScaleDisplay(this,viewer)

spacing = 25;
c = get(viewer.Figure,'Color');

p = [this.PanelXMargin + 4*this.TextBoxWidth + 2*spacing - 100,...
     this.PanelYMargin + this.TextBoxHeight-3,...
     this.TextBoxWidth,...
     this.TextBoxHeight];

uicontrol('Parent',this.LayoutPanel,'Style','text','String','Scale:',...
          'Units','Pixels','Enable','inactive',...
          'Position',p-[30 0 0 0],'HorizontalAlignment','left','BackgroundColor',c);

p = [this.PanelXMargin + 5*this.TextBoxWidth + 2*spacing - 110,...
     this.PanelYMargin + this.TextBoxHeight+2,...
     getMapUnitStrLen(this) + 20,...
     this.TextBoxHeight];

viewer.ScaleDisplay = uicontrol(this.LayoutPanel, ...
    'String',' ','Style','Edit', 'Units','Pixels', ...
    'Enable','inactive','SelectionHighlight','off',...
    'Position',p-[30 0 0 0],'Callback',@localSetScale);
                            
    function localSetScale(hSrc,event) %#ok<INUSD>
        err = false;
        str = get(hSrc,'String');
        % Remove commas
        str(str==',') = [];
        % Remove spaces
        str(isspace(str)) = [];
        if isempty(findstr(str,':'))
            num = 1;
            den = str2double(str);
            if isempty(den)
                err = true;
            end
        else
            [values, count] = sscanf(str,'%f:%f');
            if count ~= 2
                err = true;
            else
                num = values(1);
                den = values(2);
            end
        end
        if err
            oldscale = viewer.Axis.getScale;
            viewer.Axis.setScale(oldscale);
        else
            viewer.Axis.setScale(num/den);
        end
    end

end

%-----------------------------------------------------------------------

function createMapUnitsDisplay(this,viewer)

spacing = 25;
c = get(viewer.Figure,'Color');

p = [this.PanelXMargin + 4*this.TextBoxWidth + 2*spacing - 100,...
     this.PanelYMargin - 3,...
     this.TextBoxWidth,...
     this.TextBoxHeight];

uicontrol('Parent',this.LayoutPanel,'Style','text','String','Map units:',...
          'Units','Pixels','Enable','inactive',...
          'Position',p-[30 0 0 0],'HorizontalAlignment','left','BackgroundColor',c);

p = [this.PanelXMargin + 5*this.TextBoxWidth + 2*spacing - 110,...
     this.PanelYMargin,...
     getMapUnitStrLen(this) + 20,...
     this.TextBoxHeight];

mapunits = this.MapUnits(:,1)';
this.MapUnitsDisplay = uicontrol( ...
     'Parent', this.LayoutPanel, ...
     'Style', 'popupmenu', ...
     'Units', 'pixels', ...  
     'Position', p - [30 0 0 0], ...
     'String', mapunits, ... 
     'HorizontalAlignment', 'left', ...
     'Callback', @(hSrc, event) localSetMapUnits(hSrc, event, this, viewer));
end

%-----------------------------------------------------------------------

function createActiveLayerDisplay(this,viewer)

spacing = 25;
c = get(viewer.Figure,'Color');

extraPadding = this.TextBoxWidth - (getMapUnitStrLen(this) + 20);

p = [this.PanelXMargin + 6*this.TextBoxWidth + 3*spacing - 110 - extraPadding,...
     this.PanelYMargin + this.TextBoxHeight-3,...
     this.TextBoxWidth + 20,...
     this.TextBoxHeight];

uicontrol(this.LayoutPanel,'Style','text','Units','Pixels', ...
          'Position',p,'String','Active layer:',...
          'HorizontalAlignment','left','BackgroundColor',c);

p = [this.PanelXMargin + 6*this.TextBoxWidth + 3*spacing - 110 - extraPadding,...
     this.PanelYMargin,...
     getMapUnitStrLen(this) + 20,...
     this.TextBoxHeight];
%     this.TextBoxWidth + 25,... % Make it a little bigger to hold long names

this.ActiveLayerDisplay = uicontrol(this.LayoutPanel,'Style', ...
                                    'popupmenu','Units','Pixels', ...
                                    'Position',p,'String',...
                                    [{'None'};viewer.getMap.getLayerOrder],...
                                    'HorizontalAlignment','left',...
                                    'Callback',{@localSetActiveLayer this viewer});
end

%-----------------------------------------------------------------------
     
function localSetActiveLayer(hSrc,event,this,viewer) %#ok<INUSL>
val = get(this.ActiveLayerDisplay,'Value');
strs = get(this.ActiveLayerDisplay,'String');
activelayer = strs{val};
setActiveLayer(viewer,activelayer);
end

%-----------------------------------------------------------------------

function localSetMapUnits(hSrc,event,this,viewer) %#ok<INUSL>
mapUnitsInd = get(this.MapUnitsDisplay,'Value');
mapUnitsTag = this.MapUnits{mapUnitsInd,2};

mapUnitsMenu = findobj(get(viewer.Figure,'Children'),'Label','Set Map Units');
selectedMapUnitMenu = findobj(get(mapUnitsMenu,'Children'),'Tag', mapUnitsTag);
set(get(mapUnitsMenu,'Children'),'Checked','off');
set(selectedMapUnitMenu,'Checked','on');

viewer.setMapUnits(mapUnitsTag)
end

%-----------------------------------------------------------------------

function lenOut = getMapUnitStrLen(this)
% Get length of longest Map Unit
maxCharLen = length(strvcat(this.MapUnits{:,1})); 

% convert to size in pixels
lenOut = maxCharLen*5 + 10;
end
