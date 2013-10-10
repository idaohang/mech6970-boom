function newLine = makeCopy(this)
%makeCopy Copy this DragLine object
%
%   makeCopy returns a new drag line that is identical to this object,
%   except that the copy is invisible.

% Copyright 2008 The MathWorks, Inc.
% $Revision: 1.1.6.2 $  $Date: 2008/11/24 14:59:07 $

% Construct an invisible HG line object with nominal properties.
ax = get(this.hLine,'Parent');
newLine = line(...
    'XData', get(this.hLine,'XData'),...
    'YData', get(this.hLine,'YData'),...
    'Parent',ax,...
    'Selected','off',...
    'Visible','off');

% Assign the rest of the line properties.
readOnlyProperties = {'Annotation','BeingDeleted','Type'};
props = rmfield(get(this.hLine),...
    [readOnlyProperties {'Selected','Visible'}]);
set(newLine,props)

% Construct a new default MapGraphics.DragLine object.
hCopy = MapGraphics.DragLine;

% Assign selected properties.
hCopy.hLine = newLine;
hCopy.IsArrow = this.IsArrow;
hCopy.Finished = true;

% Set back-pointer from HG line object to MapGraphics.DragLine object.
setappdata(hCopy.hLine,'AnnotationObject',hCopy)

% If its HG line object is deleted, also delete this DragLine object.
set(hCopy.hLine, 'DeleteFcn', @deleteDragLine)

    function deleteDragLine(hSrc,event)  %#ok<INUSD>
        if ishandle(hCopy) % Applying ishandle to a non-HG object
            delete(hCopy)
        end
    end
end
