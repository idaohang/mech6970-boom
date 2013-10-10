function newText = makeCopy(this)
%makeCopy Copy this Text object
%
%   makeCopy constructs a new MapGraphics.Text object that is identical
%   this one, and that contains a handle to a new HG text object which
%   is identical to the original except that it is invisible and
%   unselected. The return value is handle to the new HG text object. 

% Copyright 1996-2008 The MathWorks, Inc.
% $Revision: 1.1.6.5 $  $Date: 2008/10/26 14:25:53 $

% Construct an invisible HG text object with nominal properties.
parent = get(this.hText,'Parent');
newText = text(...
    'Position', get(this.hText,'Position'),...
    'String',   get(this.hText,'String')',...
    'Parent',   parent,...
    'Selected', 'off',...
    'Visible',  'off');

% Assign the rest of the properties.
readOnlyProperties = {'Annotation','BeingDeleted','Extent','Type'};
props = rmfield(get(this.hText),...
    [readOnlyProperties {'Selected','Visible'}]);
set(newText,props)

% Construct a new default MapGraphics.Text object.
hCopy = MapGraphics.Text;

% Assign hText to it.
hCopy.hText = newText;

% Set back-pointer from HG text object to MapGraphics.Text object.
setappdata(hCopy.hText,'AnnotationObject',hCopy)

% If its HG text object is deleted, also delete this Text object.
set(hCopy.hText, 'DeleteFcn', @deleteText)

    function deleteText(hSrc,event)  %#ok<INUSD>
        if ishandle(hCopy) % Applying ishandle to a non-HG object
            delete(hCopy)
        end
    end
end
