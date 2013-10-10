function this = DefaultState(viewer)

% Copyright 1996-2008 The MathWorks, Inc.
% $Revision: 1.1.6.4 $  $Date: 2008/11/24 14:59:53 $

this = MapViewer.DefaultState;

viewer.setDefaultWindowButtonFcn();

% Menus
set(viewer.NewViewAreaMenu,'Enable','off')
set(viewer.ExportAreaMenu, 'Enable','off')
