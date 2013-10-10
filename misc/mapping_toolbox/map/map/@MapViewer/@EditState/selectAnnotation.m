function selectAnnotation(this, hObj)

%   Copyright 2005-2007 The MathWorks, Inc.
%   $Revision: 1.1.6.2 $  $Date: 2007/12/10 21:38:51 $

hAx = this.AnnotationAxes;

if isempty(hObj)
    return
end

% unselect currently selected annotations
allChildObjs = get(hAx,'Children');

unselectedObjs = setxor(hObj, allChildObjs);

nonSelectableAnnotations = findobj(allChildObjs,'tag','XLabel','-or',...
                                       'tag','YLabel','-or',...
                                       'tag','Title');

set(hObj,'Selected','on');
set(unselectedObjs,'Selected','off');
set(nonSelectableAnnotations,'Selected','off');
this.enableMenus;
