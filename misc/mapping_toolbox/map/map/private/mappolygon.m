function varargout = mappolygon(xdata, ydata, varargin)
%MAPPOLYGON Display polygon without projection
%
%   Construct a Polygon object for display in map (x-y) coordinates.
%   xdata and ydata contain the polygon vertices, and use NaNs to
%   designate multiple parts, including inner rings.  XDATA and YDATA
%   must be the same size and have NaNs in matching locations.  Polygons
%   support the same graphics properties as patches.  Return empty if
%   XDATA and YDATA are empty.
%
%   Example
%   -------
%   load coast
%   figure
%   h = mappolygon(long, lat)
%   class(h)
%   set(h,'FaceColor',[0.7 0.7 0.4])

% Copyright 2006-2009 The MathWorks, Inc.
% $Revision: 1.1.6.8 $  $Date: 2009/11/09 16:26:06 $

% Separate out any 'Parent' properties from varargin
qParent = strncmpi(varargin,'pa',2);
qParent = qParent | circshift(qParent,[0 1]);
parents = varargin(qParent);
varargin(qParent) = [];

if ~isempty(xdata) || ~isempty(ydata)
    if any(~isnan(xdata(:))) || any(~isnan(ydata(:)))
        % Clean up data, making sure that the edge-line closes.
        [xdata, ydata] = closePolygonParts(xdata, ydata);

        if isShapeMultipart(xdata,ydata)
            [f,v] = polygonToFaceVertex(xdata,ydata);
            geodata = {'Faces', f, 'Vertices', v};
        else
            % Remove  trailing NaNs.
            xdata(isnan(xdata)) = [];
            ydata(isnan(ydata)) = [];
            geodata = {'XData', xdata, 'YData', ydata};
        end
    else
        % xdata and ydata both contain nothing but NaN.
        geodata = {'XData', NaN, 'YData', NaN};
    end
    % Construct a pale-yellow patch with edges turned off;
    % keep it invisible for now.
    defaultFaceColor = [1 1 0.5];   % Pale yellow
    h = patch(geodata{:}, parents{:}, ...
        'FaceColor',defaultFaceColor, ...
        'EdgeColor','none', ...
        'Visible','off');
    
    % Construct an "edge line" object.
    hEdgeLine = constructEdgeLine(h,xdata,ydata);
    
    % After setting up listeners and initializing the "update" state, we're
    % ready to set the rest of the input properties and make both patch and
    % line visible. Also add a delete function to clean up the edge line.
    set(h,'Visible','on',varargin{:},'DeleteFcn',@deleteEdgeLine)
    set(hEdgeLine,'Visible','on')    
else
    % xdata and ydata are both empty.
    h = reshape([],[0 1]);
end

% Suppress output if called with no return value and no semicolon.
if nargout > 0
    varargout{1} = h;
end

    %------------------- nested callback function ------------------
    
    function deleteEdgeLine(hSrc,evnt) %#ok<INUSD>
        % Delete the edge line.
        if ishghandle(hEdgeLine)
            delete(hEdgeLine);
        end
    end

end
