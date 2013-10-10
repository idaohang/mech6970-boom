function hEdgeLine = constructEdgeLine(h,xdata,ydata,zdata)
% Construct an "edge line" object and set up listeners and state
% information to re-direct certain property settings.

% Copyright 2009-2010 The MathWorks, Inc.
% $Revision: 1.1.6.3 $  $Date: 2010/06/07 16:34:06 $

parent = get(h, 'Parent');
defaultColor = [0 0 0]; % Black
if nargin < 4
    hEdgeLine = line('Parent', parent, 'XData', xdata, 'YData', ydata, ...
        'Color',defaultColor,'Visible','off','HandleVisibility','off');
else
    hEdgeLine = line('Parent', parent, ...
        'XData', xdata, 'YData', ydata, 'Zdata', zdata, ...
        'Color',defaultColor,'Visible','off','HandleVisibility','off');
end

% Except for the first one (which forces the patch EdgeAlpha property
% to remain at a value of 1), each listener has a callback that re-maps a
% certain patch property to an edge line property.
addlistener(h, 'EdgeAlpha',       'PostSet', @setEdgeAlpha);
addlistener(h, 'EdgeColor',       'PostSet', @setEdgeColor);
addlistener(h, 'LineWidth',       'PostSet', @setEdgeLineProps);
addlistener(h, 'LineStyle',       'PostSet', @setEdgeLineProps);
addlistener(h, 'Visible',         'PostSet', @setEdgeLineProps);
addlistener(h, 'Parent',          'PostSet', @setEdgeLineProps);
addlistener(h, 'Marker',          'PostSet', @setMarkerProps);
addlistener(h, 'MarkerEdgeColor', 'PostSet', @setMarkerProps);
addlistener(h, 'MarkerFaceColor', 'PostSet', @setMarkerProps);
addlistener(h, 'MarkerSize',      'PostSet', @setMarkerProps);

% Keep some state information to help the listener callbacks work.
update.EdgeAlpha = true;
update.EdgeLineColor = true;
update.MarkerProps = true;

% Cache hEdgeLine in the patch's appdata to support automated tests.
setappdata(h,'EdgeLine',hEdgeLine)

    %------------------- nested callback functions ------------------
    
    function setEdgeAlpha(hSrc,evnt) %#ok<INUSD>
        % Ensure that the 'EdgeAlpha' patch property is always unity.
        if update.EdgeAlpha
            update.EdgeAlpha = false;
            set(h,'EdgeAlpha',1)
        end
        update.EdgeAlpha = true;
    end

    function setEdgeColor(hSrc,evnt) %#ok<INUSL>
        % Apply the 'EdgeColor' value to the edge line rather
        % than to the edge of every face in the patch.
        if update.EdgeLineColor
            color = evnt.AffectedObject.EdgeColor;
            if ~any(strncmpi(color,{'flat','interp'},numel(color)))
                % Filter out color values that match 'flat' or 'interp'.
                % These values can used to set the EdgeColor property of
                % a patch but not the Color property of a line.
                set(hEdgeLine, 'Color', color)
            end
            update.EdgeLineColor = false;
            set(h,'EdgeColor','none')
        end
        update.EdgeLineColor = true;
    end

    function setEdgeLineProps(hSrc,evnt)
        % Apply 'LineWidth','LineStyle', 'Visible', and 'Parent' property
        % settings to both the patch to the edge line.
        hPatch = evnt.AffectedObject;
        set(hEdgeLine, hSrc.Name, get(hPatch, hSrc.Name))
    end

    function setMarkerProps(hSrc,evnt)
        % Remap marker properties from the patch to the edge line.
        if update.MarkerProps
            hPatch = evnt.AffectedObject;
            set(hEdgeLine, hSrc.Name, get(hPatch, hSrc.Name))
            update.MarkerProps = false;
            set(h,'Marker','none')
        end
        update.MarkerProps = true;
    end
end
