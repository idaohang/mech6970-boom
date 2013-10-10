function this = InsertArrowState(viewer)

% Copyright 1996-2008 The MathWorks, Inc.
% $Revision: 1.1.6.7 $  $Date: 2008/10/26 14:26:19 $

this = MapViewer.InsertArrowState;
this.MapViewer = viewer;
viewer.setCursor({'Pointer','crosshair'});
set(viewer.Figure, 'WindowButtonDownFcn', @insertLine);

    function insertLine(hSrc,event) %#ok<INUSD>        
        if viewer.isOverMapAxes()
            MapGraphics.DragLine(viewer.AnnotationAxes,true);
        end
    end
end
