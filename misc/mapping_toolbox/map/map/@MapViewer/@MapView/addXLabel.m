function addXLabel(this)

% Copyright 1996-2008 The MathWorks, Inc.
% $Revision: 1.1.6.8 $  $Date: 2008/10/26 14:26:29 $

prompt = 'Enter the xlabel:';
dlgTitle = 'XLabel';
lineNo = 1;
if isempty(this.XLabel)
  answer=inputdlg(prompt,dlgTitle,lineNo);
  if length(answer) >= 1 && ~isempty(answer{1})
    nLines = size(answer{1},1);
    space = 14 * nLines;
    
    % Move the axis up and shrink the axis height
    this.ExtraBottomMargin = space;
    this.Axis.adjustPositionInPoints([0 space 0 -space])
        
    % X,Y text position for Units = normalized 
    xpos = 0.5;
    ypos = 0;
    
    h = text(...
        'Parent',this.AnnotationAxes,...
        'String',cellstr(answer),...
        'Color','black',...
        'HitTest','off',...
        'HorizontalAlignment','center',...
        'VerticalAlignment','top',...
        'Units','normalized',...
        'Position',[xpos ypos 0],...
        'Tag','XLabel',...
        'DeleteFcn',{@deleteXLabel this nLines});
    
    this.XLabel = h;
  end
else
  answer=inputdlg(prompt,dlgTitle,lineNo,get(this.XLabel,'String'));
  if  length(answer) >= 1 && ~isempty(answer{1})
    % Use CELLSTR to handle multi-line titles.
    set(this.XLabel,'String',cellstr(answer))
  end
end

function deleteXLabel(hSrc,event,this,nLines) %#ok<INUSL>
% Move the axis down 14 points and expand the axis height by 14 points

space = nLines * 14;
this.Axis.adjustPositionInPoints([0 -space 0 space])
this.ExtraBottomMargin = 0;
