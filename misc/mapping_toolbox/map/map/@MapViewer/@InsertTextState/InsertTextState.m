function this = InsertTextState(viewer)

% Copyright 1996-2008 The MathWorks, Inc.
% $Revision: 1.1.6.9 $  $Date: 2008/10/26 14:26:22 $

this = MapViewer.InsertTextState;
this.MapViewer = viewer;
viewer.setCursor({'Pointer','ibeam'});
set(viewer.Figure,'WindowButtonDownFcn',@doTextInsert);

    function doTextInsert(hSrc,event) %#ok<INUSD>        
        if viewer.isOverMapAxes()
            pt = get(viewer.getAxes(),'CurrentPoint');
            MapGraphics.Text(...
                'Position',[pt(1),pt(3),0],...
                'String',' ',...
                'Editing','on',...
                'Parent',viewer.AnnotationAxes);
        end
    end
end
