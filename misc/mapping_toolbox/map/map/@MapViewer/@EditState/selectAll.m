function selectAll(this)

%   Copyright 2005 The MathWorks, Inc.
%   $Revision: 1.1.6.1 $  $Date: 2005/03/31 16:33:46 $

hAx = this.AnnotationAxes;
annotations = findobj(get(hAx,'Children'),'Visible','on');

this.selectAnnotation(annotations);








