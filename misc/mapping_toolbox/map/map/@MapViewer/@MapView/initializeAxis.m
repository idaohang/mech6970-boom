function ax = initializeAxis(this,map,visibility)
%INITIALIZEAXIS Initialize axis to view the map.

%   Copyright 1996-2006 The MathWorks, Inc.
%   $Revision: 1.1.6.7 $  $Date: 2006/12/10 20:04:33 $

borderWidth = 2;
this.AxisPanel = uipanel('Parent',this.Figure,...
                         'Units','Pixels',...
                         'Position',getPanelPosition(this.Figure),...
                         'BorderType','none',...
                         'BorderWidth',borderWidth);

ax = MapGraphics.axes('Parent',this.AxisPanel,...
                      'Units','normalized','Position',[0 0 1 1],...
                      'XTick',[],'YTick',[],'XColor','white','YColor','white',...
                      'Visible',visibility);

ax.addModel(map);


function pixPanelPos = getPanelPosition(fig)
p = getpixelposition(fig);

leftMargin = 10;
bottomMargin = 65;
topMargin = 10;
rightMargin = 10;

pixPanelPos = [leftMargin,bottomMargin,p(3)-leftMargin-rightMargin,...
               p(4) - bottomMargin - topMargin];
