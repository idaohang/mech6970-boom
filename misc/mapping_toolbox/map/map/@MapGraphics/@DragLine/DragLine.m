function this = DragLine(ax,isArrow)
%DRAGLINE Draw a line on an axis between two points.

% Copyright 1996-2008 The MathWorks, Inc.
% $Revision: 1.1.6.10 $  $Date: 2008/11/24 14:59:06 $

% Construct a default DragLine object.
this = MapGraphics.DragLine;

% Construct an HG line object (with no data) and keep track of its handle.
this.hLine ...
    = line('XData',[],'YData',[],'ZData',[],'Parent',ax,'Visible','off');

% Define additional variables needed in the nested stopDrag callback.
hFig = ancestor(ax,'Figure');
oldWindowButtonMotionFcn = get(hFig,'WindowButtonMotionFcn'); 

% Initialize DragLine properties.
this.Finished = false;
this.IsArrow = isArrow;
this.OldPointer = get(hFig,'Pointer');

% Initialize DragLine starting point.
set(hFig,'CurrentObject',ax)
pt = get(ax,'CurrentPoint');
this.StartX = pt(1);
this.StartY = pt(3);

% Assign starting point coordinates to the HG line object.
set(this.hLine,'XData',this.StartX,'YData',this.StartY,'Visible','on')

% If its HG line object is deleted, also delete this DragLine object.
set(this.hLine, 'DeleteFcn', @deleteDragLine)

% Set back-pointer from the HG line object to this DragLine object.
setappdata(this.hLine,'AnnotationObject',this)

% Set up figure callbacks to support selection of end point.
set(hFig,'WindowButtonMotionFcn',@lineDrag,'WindowButtonUpFcn',@stopDrag)

%------------------ nested callback functions -------------------------

    function lineDrag(hSrc,event)  %#ok<INUSD>
        pt = get(ax,'CurrentPoint');
        set(this.hLine, ...
            'XData', [this.StartX pt(1)], ...
            'YData', [this.StartY pt(3)]);
    end

    function stopDrag(hSrc,event)  %#ok<INUSD>        
        pt = get(ax,'CurrentPoint');
        this.EndX = pt(1);
        this.EndY = pt(3);        
        set(this.hLine, ...
            'XData', [this.StartX this.EndX], ...
            'YData', [this.StartY this.EndY]);
        
        set(hFig,'CurrentObject',ax)
        
        if this.IsArrow
            this.ArrowHead = MapGraphics.ArrowHead(this.hLine);
        end
        
        set(hFig,'Pointer',this.OldPointer)
        iptPointerManager(hFig,'Enable');
        
        toolbar = findall(hFig,'type','uitoolbar');
        toolButton = findall(toolbar,'ToolTipString','Insert Line');
        
        set(toolButton,'State','off')
        this.Finished = true;
        
        % This keeps things going ... note that we are depending on the
        % function-scoped value of isArrow.  In many cases we could also
        % use this.isArrow, but the problem with that is that it
        % requires 'this' to exist at the time the button down event
        % occurs ... and reading down a few more lines you can see that
        % the DragLine handle 'this' gets deleted if a degenerate point
        % is detected.
        set(hFig,'WindowButtonDownFcn', ...
            @(hSrc,event) MapGraphics.DragLine(ax,isArrow))
        
        set(hFig,'WindowButtonMotionFcn',oldWindowButtonMotionFcn)
        set(hFig,'WindowButtonUpFcn','')
        
        % Destroy line annotations if they are zero length.
        degeneratePoint = ...
            (this.EndX == this.StartX) &&...
            (this.EndY == this.StartY);
        
        if degeneratePoint
            delete(this.hLine)
        end
    end

    function deleteDragLine(hSrc,event)  %#ok<INUSD>
        if ~isempty(this.ArrowHead) ...
                && ishandle(this.ArrowHead) ... % Applying ishandle to a non-HG object
                && ishghandle(this.ArrowHead.hPatch)
            delete(this.ArrowHead.hPatch)
        end
        if ishandle(this) % Applying ishandle to a non-HG object
            delete(this)
        end
    end
end
