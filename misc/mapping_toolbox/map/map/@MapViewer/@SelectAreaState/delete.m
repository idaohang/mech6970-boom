function delete(this)

% Copyright 1996-2008 The MathWorks, Inc.
% $Revision: 1.1.6.7 $  $Date: 2008/11/24 15:00:12 $

viewer = this.MapViewer;
set(viewer.Figure,'WindowButtonDownFcn','');
ax = viewer.getAxes();
set(ax,'ButtonDownFcn','')
set(get(ax,'Children'),'ButtonDownFcn','')

if ~isempty(this.Box) && ishghandle(this.Box)
  delete(this.Box);
end
this.disableMenus;
