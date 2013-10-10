function deleteAnnotation(this)
% Delete any currently-selected DragLines and Text objects from the
% annotation axes.

% Copyright 2005-2008 The MathWorks, Inc.
% $Revision: 1.1.6.2 $  $Date: 2008/10/26 14:26:15 $

selectedObjects = this.getSelectedAnnotation;
delete(selectedObjects)
this.enableMenus();
