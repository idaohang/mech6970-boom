function [width, height] = getVisibleSizeInPixels(this)
% Visible size of axes in pixels, allowing for text labels.

% Copyright 2008 The MathWorks, Inc.
% $Revision: 1.1.6.2 $  $Date: 2008/11/24 14:59:28 $

ax = this.getAxes();

% Lower Right
t = text('Parent',ax,'Units','normalized','Position',[0 0 0],...
         'String','ignore','Visible','off');
set(t,'Units','pixels');
p = get(t,'Position');
px1 = p(1);
py1 = p(2);
delete(t);

% Upper Left
t = text('Parent',ax,'Units','normalized','Position',[1 1 0],...
         'String','ignore','Visible','off');
set(t,'Units','pixels');
p = get(t,'Position');
px2 = p(1);
py2 = p(2);
delete(t);

width  = px2 - px1;
height = py2 - py1;
