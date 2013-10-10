function hLine = cut(this)
%CUT Cut this DragLine object

% Copyright 2008 The MathWorks, Inc.
% $Revision: 1.1.6.1 $  $Date: 2008/10/26 14:25:41 $

% Note: If the drag line has an ArrowHead object, its setVisible
%       listener will the visibility of its patch in synch with that of
%       the HG line.

hLine = this.hLine;
set(hLine,'Selected','off','Visible','off')
