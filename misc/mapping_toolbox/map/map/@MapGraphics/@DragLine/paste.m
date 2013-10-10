function paste(this,xShift,yShift)

% Copyright 2008 The MathWorks, Inc.
% $Revision: 1.1.6.2 $  $Date: 2008/11/24 14:59:08 $

xData = get(this.hLine,'XData');
yData = get(this.hLine,'YData');

set(this.hLine,...
    'XData',xData + xShift,...
    'YData',yData + yShift,...
    'Visible','on')

% If this DragLine has an ArrowHead, then it will be visible now because
% its patch component listens for changes in the visibility of
% this.hLine.  This will be the case if we're pasting a DragLine (with
% arrow head) that was cut. But if we're pasting a copy, it doesn't have
% an arrow head yet (because we've waited to create it until we need it,
% which is right now).
if this.IsArrow
    if isempty(this.ArrowHead) ...
            || ~ishandle(this.ArrowHead) ... % Applying ishandle to a non-HG object
            || ~ishghandle(this.ArrowHead.hPatch,'patch')
        this.ArrowHead = MapGraphics.ArrowHead(this.hLine);
    end
end
