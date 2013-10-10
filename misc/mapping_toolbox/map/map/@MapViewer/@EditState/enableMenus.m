function enableMenus(this)

% Copyright 1996-2008 The MathWorks, Inc.
% $Revision: 1.1.6.6 $  $Date: 2008/10/26 14:26:16 $

% set(get(this.EditMenu,'Children'),'Enable','on');

cutMenu = findall(this.EditMenu,'Tag','cut');
copyMenu = findall(this.EditMenu,'Tag','copy');
pasteMenu = findall(this.EditMenu,'Tag','paste');
selectAllMenu = findall(this.EditMenu,'Tag','select all');

viewer = this.MapViewer;

allChildObjs = get(viewer.AnnotationAxes,'Children');

nonSelectableAnnotations = findobj(allChildObjs,'tag','XLabel','-or',...
                                       'tag','YLabel','-or',...
                                       'tag','Title');

isAnnotation = ~isempty(allChildObjs) && ~isempty(setdiff(allChildObjs,nonSelectableAnnotations));

isCopiedObject = any(cellfun(...
    @ishghandle,viewer.CopiedObjects));

isObjectSelected = ~isempty(this.getSelectedAnnotation);

% handle any annotation displayed 
if isAnnotation
    set(selectAllMenu,'Enable','on');
else
    set(selectAllMenu,'Enable','off');
end

% handle copy and cut menu items
if isObjectSelected 
    set([cutMenu,copyMenu],'Enable','on');
else
    set([cutMenu,copyMenu],'Enable','off');
end

% handle paste menu
if isCopiedObject
    set(pasteMenu,'Enable','on');
else
    set(pasteMenu,'Enable','off');
end
