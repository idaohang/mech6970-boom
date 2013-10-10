function unselectAnnotation(this, hObj)

%   Copyright 2005 The MathWorks, Inc.
%   $Revision: 1.1.6.1 $  $Date: 2005/03/31 16:33:48 $

hAx = this.AnnotationAxes;
if nargin == 1
    annotations = get(hAx,'Children');
else
    annotations = hObj;
end
set(annotations,'Selected','off');
this.enableMenus;
