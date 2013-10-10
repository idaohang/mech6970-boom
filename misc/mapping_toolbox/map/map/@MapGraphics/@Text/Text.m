function this = Text(varargin)
%TEXT 

% Copyright 1996-2008 The MathWorks, Inc.
% $Revision: 1.1.6.3 $  $Date: 2008/10/26 14:25:50 $

% Construct default Text object.
this = MapGraphics.Text;

% Construct an HG text object and keep track of its handle.
this.hText = text(varargin{:});

% Set back-pointer from the HG text object to this Text object.
setappdata(this.hText,'AnnotationObject',this)

% If its HG text object is deleted, also delete this Text object.
set(this.hText, 'DeleteFcn', @deleteText)

    function deleteText(hSrc,event)  %#ok<INUSD>
        if ishandle(this) % Applying ishandle to a non-HG object
            delete(this)
        end
    end
end
