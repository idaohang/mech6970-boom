function addTitle(this)

% Copyright 1996-2008 The MathWorks, Inc.
% $Revision: 1.1.6.8 $  $Date: 2008/10/26 14:26:28 $

prompt = 'Enter the title:';
dlgTitle = 'Title';
lineNo = 1;
if isempty(this.Title)
  answer=inputdlg(prompt,dlgTitle,lineNo);
  if length(answer) >= 1 && ~isempty(answer{1})
    nLines = size(answer{1},1);
    space = 12 * nLines;
    this.ExtraTopMargin = space;
    this.Axis.adjustPositionInPoints([0 0 0 -space])
    
    % X,Y text position for Units = normalized 
    xpos = 0.5;
    ypos = 1.0;
    
    h = text(...
        'Parent',this.AnnotationAxes,...
        'String',cellstr(answer),...
        'Color','black',...
        'HitTest','off',...
        'HorizontalAlignment','center',...
        'VerticalAlignment','bottom',...
        'Units','normalized',...
        'Position',[xpos ypos 0],...
        'Tag','Title',...
        'DeleteFcn',{@deleteTitle this nLines});
    
    this.Title = h;
  end
else
  answer=inputdlg(prompt,dlgTitle,lineNo,get(this.Title,'String'));
  if  length(answer) >= 1 && ~isempty(answer{1})
    % Use CELLSTR to handle multi-line titles.
    set(this.Title,'String',cellstr(answer))
  end
end

function deleteTitle(hSrc,event,this,nLines) %#ok<INUSL>
% Move the axis down 12 points and expand the axis height by 12 points

space = nLines * 12;
this.Axis.adjustPositionInPoints([0 0 0 space]);
this.Title = [];
this.ExtraTopMargin = 0;
