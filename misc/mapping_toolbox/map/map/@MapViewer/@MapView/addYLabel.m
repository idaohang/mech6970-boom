function addYLabel(this)

% Copyright 1996-2008 The MathWorks, Inc.
% $Revision: 1.1.6.7 $  $Date: 2008/10/26 14:26:30 $

if isempty(this.YLabel)
  answer=inputdlg('Enter the ylabel:','YLabel',1);
  if  length(answer) >= 1 && ~isempty(answer{1})
    nLines = size(answer{1},1);
    space = 12 * nLines;
    % Move the axis right and shrink the axis width
    this.ExtraLeftMargin = space;
    this.Axis.adjustPositionInPoints([space 0 -space 0])
    
    % X,Y text position for Units = normalized 
    xpos = 0;
    ypos = 0.5;
    
    h = text(...
        'Parent',this.AnnotationAxes,...
        'String',cellstr(answer),...
        'Color','black',...
        'HitTest','off',...
        'HorizontalAlignment','center',...
        'VerticalAlignment','bottom',...
        'Units','normalized',...
        'Position',[xpos ypos 0],...
        'Rotation',90,...
        'Tag','YLabel',...
        'DeleteFcn',{@deleteYLabel this nLines});
    
    this.YLabel = h;
  end
else
  answer=inputdlg('Enter the YLabel:','YLabel',2,get(this.YLabel,'String'));
  if  length(answer) >= 1 && ~isempty(answer{1})
    % Use CELLSTR to handle multi-line titles.
    set(this.YLabel,'String',cellstr(answer))
  end
end

function deleteYLabel(hSrc,event,this,nLines) %#ok<INUSL>
% Move the axis over 12 points and expand the axis width by 12 points

space = nLines * 12;
this.Axis.adjustPositionInPoints([-space 0 space 0])
this.ExtraLeftMargin = 0;
