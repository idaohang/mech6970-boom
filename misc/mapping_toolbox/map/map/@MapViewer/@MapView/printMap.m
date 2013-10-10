function printMap(this)

% Copyright 1996-2008 The MathWorks, Inc.
% $Revision: 1.1.6.4 $  $Date: 2008/11/24 15:00:06 $

fig = figure('Visible','off','Position',get(this.Figure,'Position'));
copyobj(this.getAxes(),fig);
% Avoid deleting this.Axis (a MapGraphics.axes object) when closing fig
% -- removed the delete function from the axes we've just cloned.
hAx = findobj(fig,'Type','axes');
set(hAx,'DeleteFcn',[])
set(this.AnnotationAxes,'Parent',fig)
axis 'off';
printdlg(fig);
set(this.AnnotationAxes,'parent',this.Figure)
close(fig);
