function selObj = getSelectedAnnotation(this)

%   Copyright 2005 The MathWorks, Inc.
%   $Revision: 1.1.6.1 $  $Date: 2005/03/31 16:33:44 $

hAx = this.AnnotationAxes;
selObj = findobj(get(hAx,'Children'),'Selected','on');

