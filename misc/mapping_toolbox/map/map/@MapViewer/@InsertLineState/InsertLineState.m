function this = InsertLineState(viewer)

% Copyright 1996-2008 The MathWorks, Inc.
% $Revision: 1.1.6.8 $  $Date: 2008/10/26 14:26:20 $

this = MapViewer.InsertLineState;
this.MapViewer = viewer;
viewer.setCursor({'Pointer','crosshair'});
set(viewer.Figure,'WindowButtonDownFcn',@insertLine);

    function insertLine(hSrc,event) %#ok<INUSD>
        if viewer.isOverMapAxes()
            MapGraphics.DragLine(viewer.AnnotationAxes,false);
        end
    end
end
