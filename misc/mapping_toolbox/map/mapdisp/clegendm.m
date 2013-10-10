function hLegend = clegendm(varargin)
%CLEGENDM Add legend labels to map contour display
%
%   CLEGENDM(CS, H) adds a legend specifying the contour line heights to the
%   current map contour plot.  CS and H are the contour matrix output and
%   object handle outputs from CONTOURM, CONTOUR3M, or CONTOURFM.
%
%   CLEGENDM(CS, H, LOC) places the legend in the specified location:
%
%        0 = Automatic placement (default)
%        1 = Upper right-hand corner
%        2 = Upper left-hand corner
%        3 = Lower left-hand corner
%        4 = Lower right-hand corner
%       -1 = To the right of the plot
%
%   CLEGENDM(...,UNITSTR) appends the character string UNITSTR to each entry
%   in the legend.
%
%   CLEGENDM(...,STRINGS) uses the strings specified in cell array STRINGS.
%   STRINGS must have same number of entries as the line children of H.
%
%   H = CLEGENDM(...) returns the handle to the LEGEND object created.
%
%   Examples
%   --------
%   % Create a legend in the upper left-hand corner with a unit string
%   % indicating that the contour elevations are in meters.
%   load topo
%   figure('Color','w'); axesm robinson; framem; tightmap; axis off
%   [cs,h] = contourm(topo,topolegend,3);
%   clegendm(cs, h, 2, ' m')
%
%   % Create a legend in the upper left-hand corner with specified strings.
%   load topo
%   figure('Color','w'); axesm robinson; framem;  tightmap; axis off
%   [cs,h] = contourm(topo,topolegend,3);
%   str = {'low altitude','medium altitude','high altitude'};
%   clegendm(cs,h,2,str)
%
%   See also CLABELM, CONTOURFM, CONTOURM, CONTOUR3M, LEGEND.

% Copyright 1996-2010 The MathWorks, Inc.
% $Revision: 1.12.4.11 $  $Date: 2010/09/13 16:16:17 $

% Validate number of inputs.
checknargin(2,4,nargin,mfilename);

% Parse the inputs.
[h, loc, strings] = parseInputs(varargin);

% Save current axes handle.  Legend creates a new axes.
axishndl = gca;

% Create the legend.
legendhndl = createLegend(h, loc, strings);

% Assign the SCRIBE callbacks to the legend objects        
children = get(legendhndl,'Children');     
set(children,'ButtonDownFcn','uimaptbx');

% Reset current axes.
set(get(axishndl,'Parent'),'CurrentAxes',axishndl)  

% Assign output if requested.
if nargout > 0
    hLegend = legendhndl;
end

%--------------------------------------------------------------------------

function [h, loc, strings] = parseInputs(inputs)
% Parse the INPUTS cell array.

% Obtain G, LOC, and STRINGS from INPUTS. Note that the contour matrix, C,
% is not needed.
switch numel(inputs)
    case 2
        % CLEGENDM(CS, H)
        % c = varargin{1};
        g = inputs{2};
        loc = 0;
        strings = {};
    case  3
        % c = varargin{1};
        g = inputs{2};
        if ischar(inputs{3}) || iscell(inputs{3})
            % CLEGENDM(CS, H, 'UNITSTR'/STRINGS)
            loc = 0;
            strings = inputs{3};
        else
            % CLEGENDM(CS, H, LOC)
           loc = inputs{3};
           strings = {};
        end
    case 4
        % CLEGENDM(CS, H, 'UNITSTR'/STRINGS)
        % c = varargin{1};
        g = inputs{2};
        loc = inputs{3};
        strings = inputs{4};
end

% Validate that the handle input is from CONTOURM, CONTOURFM, or CONTOUR3M.
assert(~isempty(g) && ishghandle(g, 'hggroup'), ...
    'map:clegendm:notHgGroupHandle', ...
    'The parameter, H, is not a hggroup handle from %s, %s, or %s.', ...
    'CONTOURM', 'CONTOURFM', 'CONTOUR3M');

% Get the GeoContourGroup object.
h = getappdata(g, 'mapgraph');

% Validate the handle.
assert(isa(h,'internal.mapgraph.GeoContourGroup'), ...
    'map:clegendm:notGeoContourGroupHandle', ...
    'The handle, H, is not a handle from %s, %s, or %s.', ...
    'CONTOURM', 'CONTOURFM', 'CONTOUR3M');

% Validate LOC
validateattributes(loc, {'numeric'}, ...
    {'scalar', 'integer', '<=',4, '>=',-1}, mfilename, 'LOC', 3);

% Convert LOC to string value
locTable = { ...
    'NorthEastOutside', ...
    'Best', ...
    'NorthEast', ...
    'NorthWest', ...
    'SouthWest', ...
    'SouthEast'};
loc = lower(locTable{loc+2});

% Validate UNITSTR, STRINGS
if ischar(strings)
    validateattributes(strings, {'char'}, ...
        {'nonempty'}, mfilename, 'UNITSTR', 4);
else
    validateattributes(strings, {'cell'}, ...
       {'2d'}, mfilename, 'STRINGS', 4);
end

%--------------------------------------------------------------------------

function legendhndl = createLegend(h, loc, strings)
% Create a legend and return its handle.

% Construct an invisible line corresponding to each contour line,
% pass the line handles to legend, then delete them.
c = get(h,'Children');
if ~isempty(c)
    contourLines = findobj(c,'Type','line');

    if ~iscell(strings)
        levelNames = cell(1,numel(contourLines));
        unitstr = strings;
    elseif ~isempty(strings)
        levelNames = strings;
    else
        levelNames = cell(1,numel(contourLines));
        unitstr = '';
    end
    
    getLevelsFromLines = ~(iscell(strings) && ~isempty(strings));
    
    for k = numel(contourLines):-1:1
        color = get(contourLines(k),'Color');
        hLine(k,1) = line([0 1], [0 1], 'Color', color, 'Visible', 'off');
        if getLevelsFromLines
            levelNames{1,k} ...
                = [num2str(getappdata(contourLines(k),'Level')) unitstr];
        end
    end
    if getLevelsFromLines
        levelNames = fliplr(levelNames);
        hLine = flipud(hLine);
    end
    legendhndl = legend(hLine, levelNames, 'Location', loc);
    delete(hLine);
else
    legendhndl = [];
end
