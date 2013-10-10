function hndlout = contourcmap(varargin)
%CONTOURCMAP Contour colormap and colorbar for current axes
% 
%   CONTOURCMAP(CMAPSTR) updates the figure's colormap for the current axes
%   with the colormap specified by the string, CMAPSTR. Valid entries for
%   CMAPSTR include  'pink', 'hsv', 'jet', or any similar MATLAB colormap
%   function. If the axes contains Mapping Toolbox contour objects, the
%   resultant colormap contains the same number of colors as the original
%   colormap. Otherwise, the resultant colormap contains 10 colors.
% 
%   CONTOURCMAP(CMAPSTR, CDELTA) updates the figure's colormap with colors
%   varying according to CDELTA. If CDELTA is a scalar, it represents a
%   step size, and colors are generated at multiples of CDELTA. If CDELTA
%   is a vector of evenly spaced values, colors are generated at those
%   values; otherwise an error is issued. If the axes contains Mapping
%   Toolbox contour objects, the value of CDELTA is ignored.
%
%   CONTOURCMAP(..., PARAM1, VAL1, PARAM2, VAL2,...) allows you to add a
%   colorbar and control the colorbar's properties. Parameter names can be
%   abbreviated and are case-insensitive. See the table below for a list of
%   available parameters.
%
%   Optional Parameters
%   -------------------
%  
%   'Colorbar'          String with values 'on' or 'off' specifying 
%                       whether a colorbar is present, 'on', or absent
%                       from the axes, 'off'.
%
%   'Location'          String specifying the location of the colorbar.
%                       Permissible values are 'vertical' (the default), 
%                       'horizontal', or 'none'.
%
%   'ColorAlignment'    String specifying the alignment of the labels in
%                       the colorbar. Permissible values are 'center',
%                       where the labels are centered on the color bands or
%                       'ends' where labels are centered on the color
%                       breaks. If the axes contains Mapping Toolbox
%                       contour objects, the ColorAlignment will be set
%                       automatically to 'center' for contour lines and
%                       'ends' for filled contours, and cannot be modified.
%
%   'SourceObject'      Handle of the graphics object which is used to 
%                       determine the color limits for the colormap. The
%                       SourceObject value is the handle of a currently
%                       displayed object. If omitted, gca is used.
%
%   'TitleString'       String specifying the title of the colorbar axes.
%
%   'XLabelString'      String specifying the X label of the colorbar axes.
%
%   'YLabelString'      String specifying the Y label of the colorbar axes.
%
%   'ZLabelString'      String specifying the Z label of the colorbar axes.
%                        
%   In addition, properties and values that can be applied to the title and
%   labels of the colorbar axes are valid.
% 
%   h = CONTOURCMAP(...) returns a handle to the colorbar axes.
% 
%   Examples
%   --------
%   load topo
%   load coast
%   figure
%   worldmap(topo, topolegend)
%   contourfm(topo, topolegend);
%   contourcmap('jet', 'Colorbar', 'on', ...
%      'Location', 'horizontal', ...
%      'TitleString', 'Contour Intervals in Meters');
%   plotm(lat, long, 'k')
%
%   load topo
%   load coast
%   figure
%   worldmap(topo, topolegend)
%   geoshow(topo, topolegend, 'DisplayType', 'texturemap');
%   contourcmap('summer', 2000, 'Colorbar', 'on', ...
%      'Location', 'horizontal', ...
%      'TitleString', 'Contour Intervals in Meters');
%   plotm(lat, long, 'k')
%
%   See also CLABELM, CLEGENDM, COLORMAP, CONTOURFM, CONTOURM.

% Copyright 1996-2010 The MathWorks, Inc.
% $Revision: 1.6.4.9 $  $Date: 2010/07/19 12:53:59 $

% Parse required parameters.
[cdelta, cmapstr, params] = parseParams(varargin{:});

% Check for parameter/value pairs.
[h, isContourHandle, cbarProps] = parseOptionalParams(params);

% Validate the special contour parameters.
[cdelta, cbarProps.coloralignment] = validateContourParams( ...
    h, isContourHandle, cdelta, cbarProps.coloralignment);

% Create a contour colormap of the object in H and return the colors and
% levels required for the colorbar.
[clevels, cindex, numcolors] = createContourColormap( ...
    h, isContourHandle, cdelta, cmapstr, cbarProps.coloralignment);

% Create or delete the colorbar.
if isequal(cbarProps.colorbar, 'on')
    % Create the colorbar.
    hCbarAxes = createColorbar(clevels, cindex, numcolors, cbarProps); 
else
    % Delete the colorbar, if found.
    deleteColorbar();
    hCbarAxes = [];
end

% Return output arguments.
if nargout > 0;
    hndlout = hCbarAxes;
end

end

%--------------------------------------------------------------------------

function [cdelta, cmapstr, params] = parseParams(varargin)

% Validate number of inputs.
checknargin(1,inf,nargin,mfilename);

% CMAPSTR
if isnumeric(varargin{1})
    checknargin(2,inf,nargin,mfilename);
    cmapstr = varargin{2};
    pos = [2,1];
    varargin(2) = [];
else
    cmapstr = varargin{1};
    varargin(1) = [];
    pos = [1,2];
end

% Validate CMAPSTR
checkinput(cmapstr, 'char', {'nonempty'}, mfilename, 'CMAPSTR', pos(1));
assert(exist(cmapstr,'file') == 2, ...
    'map:contourcmap:invalidParam', ...
    '%s must contain the name of a colormap function.','CMAPSTR')

% CDELTA
if ~isempty(varargin) && isnumeric(varargin{1})
    cdelta = varargin{1};
    checkinput(cdelta, 'numeric', {'real','nonempty'}, ...
        mfilename, 'CDELTA', pos(2));
    params = varargin(2:end);
else
    cdelta = [];
    params = varargin;
end

end

%--------------------------------------------------------------------------

function [h, isContourHandle, cbarProps] = parseOptionalParams(params)
% Parse the parameter/value pairs.

% Validate the optional parameter/value pairs.
cbarProps = validateOptionalParams(params);

% Set default colorbar location.
if isequal(cbarProps.colorbar, 'off')
    cbarProps.location = 'none';
end

% Find the source handle and determine if it is a GeoContourGroupHandle.
[h, isContourHandle] = findSourceObject(cbarProps.sourceobject);

% Assign a value to the coloralignment property, if not set.
cbarProps.coloralignment = assignDefaultColorAlignment( ...
    h, isContourHandle, cbarProps.coloralignment);

end

%--------------------------------------------------------------------------

function cbarProps = validateOptionalParams(params)
% Validate the optional parameter/value pairs.

assert(mod(numel(params),2) == 0, ...
    'map:contourcmap:invalidParams', ...
    'The property/value inputs must be supplied as pairs.');

% Assign default values.
cbarProps = struct( ...
    'titlestring', '', ...
    'xlabelstring', '', ...
    'ylabelstring', '', ...
    'zlabelstring', '', ...
    'titleParams', '', ...
    'location', 'vertical', ...
    'coloralignment', '', ...
    'colorbar', 'off', ...
    'sourceobject', gca);

% Assign values.
if ~isempty(params)
    todelete = false(size(params));
    for i=1:2:length(params)
        params{i} = canonicalProps(params{i});
        switch params{i}
            case 'location'
                cbarProps.(params{i}) = validatestring(params{i+1},  ...
                    {'none','vertical','horizontal'}, ...
                    mfilename, '''Location''');
                todelete([i i+1]) = true;
                
            case 'coloralignment'
                cbarProps.(params{i}) = validatestring(params{i+1},...
                    {'center','ends'}, mfilename, '''ColorAlignment''');
                todelete([i i+1]) = true;
                
            case 'colorbar'
                cbarProps.(params{i}) = validatestring(params{i+1}, ...
                    {'on','off'}, mfilename, '''Colorbar''');
                todelete([i i+1]) = true;
                
            case {'xlabelstring', 'ylabelstring', ...
                  'zlabelstring', 'titlestring', ...
                  'sourceobject'}
                cbarProps.(params{i}) = params{i+1};
                todelete([i i+1]) = true;
        end
    end
    params(todelete) = [];
end
cbarProps.titleParams = params;

end

%--------------------------------------------------------------------------

function [h, isContourHandle] = findSourceObject(h)
% Find the source handle from H and determine if the source handle is a
% GeoContourGroup handle.

assert(ishghandle(h), ...
    'map:contourcmap:invalidHandle', ...
    '%s must be a handle to an axes or graphic object.','''SourceObject''');

g = findobj(h, 'Type', 'hggroup');
if ~isempty(g)
    for k=1:numel(g)
        if isappdata(g(k), 'mapgraph')
            mapgraph = getappdata(g(k), 'mapgraph');
            if isa(mapgraph, 'internal.mapgraph.GeoContourGroup')
                h = mapgraph;
                isContourHandle = true;
                return
            end
        end
    end
end

% Handle type was not found.
isContourHandle = false;

end

%--------------------------------------------------------------------------

function coloralignment = assignDefaultColorAlignment(...
    h, isContourHandle, coloralignment)
% Assign the default value of coloralignment if not set and validate the
% property if it has been set.

if isempty(coloralignment)
    if isContourHandle
        if isequal(h.Fill,'on')
            coloralignment = 'ends';
        else
            coloralignment = 'center';
        end
    else
        coloralignment = 'ends';
    end
end

end

%--------------------------------------------------------------------------

function [cdelta, coloralignment] = validateContourParams( ...
    h, isContourHandle, cdelta, coloralignment)

if isContourHandle
    % Display a warning if CDELTA has been set.
    if ~isempty(cdelta)
        warning('map:contourcmap:ignoringCDELTA', ...
            'The ''%s'' parameter is ignored with contour objects.', ...
            'CDELTA');
    end
    
    % Display a warning if coloralignment is not set properly.
    if isequal(h.Fill, 'on') && ~isequal(coloralignment, 'ends')
        warning('map:contourcmap:ignoringCenterColorAlignment', ....
            ['When the contour levels are filled, the ''%s'' ', ...
            'property must be set to ''%s''.'], 'ColorAlignment', 'ends');
        coloralignment = 'ends';
        
    elseif isequal(h.Fill, 'off') && ~isequal(coloralignment, 'center')
        warning('map:contourcmap:ignoringEndsColorAlignment', ...
            ['When the contour levels are not filled, the ''%s'' ', ...
            'property must be set to ''%s''.'], 'ColorAlignment', 'center');
        coloralignment = 'center';
    end 
end

end

%--------------------------------------------------------------------------

function [clevels, cindex, numcolors] = createContourColormap( ...
    h, isContourHandle, cdelta, cmapstr, coloralignment)
% Create a contour colormap of H. 

if isContourHandle    
    % Update the colormap of the Figure, using the colormap defined by the
    % string, CMAPSTR. 
    updateContourColormap(h, isContourHandle, cmapstr)
    
    % Obtain the color levels and color index values.
    [clevels, cindex] = computeContourColors(h);
    
    % Assign the number of colors.
    numcolors = numel(cindex);
else    
    % Validate the color limits.
    climits = validateColorLimits(h);
    
    % Compute the color levels for labels with nice increments.
    clevels = computeColorLevels(cdelta, climits);
    
    % Update the colormap based on the number of colors.
    updateColormap(cmapstr, coloralignment, clevels);
    
    % Obtain the color index and the number of colors.
    cindex = 1:size(colormap,1);
    numcolors = size(cindex,1);
end

end

%--------------------------------------------------------------------------

function climits = validateColorLimits(h)
% Check that h is a handle and get color limits

if ishghandle(h,'axes')
    caxis('auto')
    climits = get(h, 'CLim');
else
    try
        cdata = get(h,'CData');
    catch e
        error('map:contourcmap:invalidHandleObject', ...
            'Object must be an axes or contain color data.')
    end
    
    if isempty(cdata) || all(isnan(cdata))
        error('map:contourcmap:invalidCDATA', ...
            'Object''s %s property must contain numeric color data.', ...
            'CData')
    end
    
    climits = [min(cdata(:)) max(cdata(:))];
    if diff(climits) == 0;
        error('map:contourcmap:invalidCDATAProperty', ...
            'Object''s %s property must contain a range of values','CData')
    end
end
  
end

%--------------------------------------------------------------------------

function clevels = computeColorLevels(cdelta, climits)
% Compute the color levels.

if isscalar(cdelta)
	cmin = floor(climits(1));
	multfactor = floor(cmin/cdelta);
	cmin = cdelta*multfactor;
	cmax = ceil(climits(2));
	clevels = cmin:cdelta:cmax;
    if max(clevels) < cmax
        clevels = [clevels max(clevels)+cdelta];
    end
elseif isempty(cdelta)
    cmin = floor(climits(1));
    cmax = ceil(climits(2));
    cdelta = abs(diff(climits))/10;
    clevels = cmin:cdelta:cmax;
else
	% check to see that spacing is the same
	tolerance = eps;
    if abs(max(diff(diff(cdelta)))) > tolerance
        error('map:contourcmap:invalidCDATASize', ...
            '%s must consist of evenly spaced elements.','CDELTA')
    end
	clevels = cdelta;
end	

% round numbers less than epsilon to zero
clevels(abs(clevels) < eps) = 0;

% Reset the climits based on new increments
set(gca, 'CLim', [min(clevels) max(clevels)])

end

%--------------------------------------------------------------------------

function updateColormap(cmapstr, coloralignment, clevels)
% Update the colormap based on the color alignment.

if isequal(coloralignment, 'ends')
    % 'ends'
    numberofcolors = length(clevels)-1;
else
    % 'center'
    numberofcolors = length(clevels);  
end
cmapfunc = str2func(cmapstr);
cmap = cmapfunc(numberofcolors);
colormap(cmap);

end

%--------------------------------------------------------------------------

function updateContourColormap(h, isContourHandle, cmapstr)
% Update the contour colormap. Keep the same number of colors as in the
% current Figure's colormap.

orig_cmap = colormap;
numberOfColors = size(orig_cmap,1);
cmapfunc = str2func(cmapstr);
cmap = cmapfunc(numberOfColors);

% Update the colormap if it has changed.
if isContourHandle
    % Update contour object
    if strcmp(h.Fill,'off')
        if ~isequal(cmap, get(h,'LineColormap'))
            set(h,'LineColormap',cmap)
            h.refresh()
        end
    else
        if ~isequal(cmap, get(h,'FillColormap'))
            set(h,'FillColormap',cmap)
            h.refresh()
        end
    end
    
    % Update figure
    f = ancestor(h.HGGroup, 'figure');
    set(f,'Colormap',cmap)
end

end

%--------------------------------------------------------------------------

function [clevels, cindex] = computeContourColors(h)
% Compute contour levels CLEVELS and a row index into the color map,
% CINDEX. For contour lines only, CINDEX includes an element for each
% element in CLEVELS. In the case of filled contours, CLEVELS includes
% an extra element equal to the upper limit of the highest contour
% interval.

if strcmp(h.Fill,'off')
    % Levels and colors for contour lines.
    S = h.getContourLines();
    levels = [S.Level];
    clevels = levels;
    cindex = internal.mapgraph.selectColors(get(h,'LineColormap'), levels);
else
    % Limits of contour intervals and colors for filling each interval.
    S = h.getFillPolygons();
    levels = [S.MinLevel];
    clevels = unique([S.MinLevel S.MaxLevel]);
    cindex = internal.mapgraph.selectColors(get(h,'FillColormap'), levels);
end

end

%--------------------------------------------------------------------------

function hAxes = createColorbar(clevels, colors, numcolors, cbarProps)
% Create the colorbar.

% Default is to use current axes.
hax = gca;
hfig = ancestor(hax,'figure');

% Get the axes units and change them to normalized.
axesunits = get(hax,'Units');
set(hax, 'Units', 'normalized')

% Create the colorbar, if present.
deleteColorbar();

% Create the structure for callbacks.
axesinfo.h = hax;
axesinfo.units = axesunits;
axesinfo.origpos = get(hax, 'Position');

switch cbarProps.location 
    case 'none'
        hAxes = [];
        
    case 'vertical'
        hAxes = createVerticalColorbar(axesinfo, hax, ...
            cbarProps.coloralignment, clevels, colors, numcolors);
        
    case 'horizontal'
        hAxes = createHorizontalColorbar(axesinfo, hax, ...
            cbarProps.coloralignment, clevels, colors, numcolors);
end

% Reset the axes units.
set(hax, 'Units', axesunits)

% Activate the initial axes.
set(hfig,'CurrentAxes',hax)

% Set text properties of colorbar axes and title if provided.
if ~isempty(hAxes)
    setColorbarProps(hAxes, cbarProps);
end

end

%--------------------------------------------------------------------------

function  hAxes = createVerticalColorbar( ...
    ax, hax, coloralignment, clevels, colors, numcolors)

xlimits = [0 1];
ylimits = [1 numcolors];

% Shrink length by 10 percent.
pos = ax.origpos;
pos(3) = ax.origpos(3)*0.90;
set(hax, 'Position', pos)

% Calculate the position of the colorbar axes.
len   = ax.origpos(3)*0.05;
width = ax.origpos(4);
axesPos = [ax.origpos(1)+ax.origpos(3)*0.95 ax.origpos(2) len width];

% Compute the tick locations.
ytickloc = getTickLocation(coloralignment, numel(clevels));

% Create the axes properties.
axesProps = struct( ...
    'Position', axesPos, ...
    'YTick', ytickloc, 'YTickLabel', clevels, ...
    'XTick', [], 'YDir', 'normal', 'YAxisLocation', 'right');

% Create the colorbar axes and image.
hAxes = colorbarAxes(axesProps, ax, xlimits, ylimits, colors');

end

%--------------------------------------------------------------------------

function hAxes = createHorizontalColorbar( ...
    ax, hax, coloralignment, clevels, colors, numcolors)

xlimits = [1,numcolors];
ylimits = [0 1];

% Shrink width by 10 percent.
pos = ax.origpos;
pos(4) = ax.origpos(4)*0.90;
pos(2) = ax.origpos(2) + ax.origpos(4)*0.10;
set(hax, 'Position', pos)

% Calculate the position of the colorbar axes.
width = ax.origpos(4)*0.05;
len = ax.origpos(3);
axesPos = [ax.origpos(1) ax.origpos(2) len width];

% Compute the tick locations.
xtickloc = getTickLocation(coloralignment, numel(clevels));

% Set the axes properties.
axesProps = struct( ...
    'Position', axesPos, ...
    'XTick', xtickloc, 'XTickLabel', clevels, ...
    'YTick', [], 'XDir', 'normal', 'XAxisLocation', 'bottom');

% Create the colorbar axes and image.
hAxes = colorbarAxes(axesProps, ax, xlimits, ylimits, colors);

end

%--------------------------------------------------------------------------

function tickloc = getTickLocation(coloralignment, numticks)
% Calculate the locations for the tick marks.

if isequal(coloralignment, 'ends')
    numticks = numticks - 1;
end

switch coloralignment
    case 'ends'
        lowerlim = 1-0.5;
        upperlim = numticks+0.5;
    case 'center'
        lowerlim = 1;
        upperlim = numticks;
end
delta = 1;
tickloc = lowerlim:delta:upperlim;

end

%--------------------------------------------------------------------------

function hAxes = colorbarAxes(axesProps, ax, xLimits, yLimits, colors)

% Create the colorbar axes.
hAxes = axes('Position', axesProps.Position);

% Create the colorbar image.
image(xLimits, yLimits, colors, ...
    'Parent', hAxes, ...
    'Tag', 'CONTOURCMAP', ...
    'DeleteFcn', @deleteim);

% Set the other axes properties.
axesProps = rmfield(axesProps, 'Position');
set(hAxes, axesProps);

% Set delete function.
setAxesDeleteFcn(hAxes, ax)

end

%--------------------------------------------------------------------------

function setColorbarProps(hndl, cbarProps)
% Set the special properties of the colorbar axes.

set(get(hndl,'Title'), 'String',cbarProps.titlestring);
set(get(hndl,'Xlabel'),'String',cbarProps.xlabelstring);
set(get(hndl,'Ylabel'),'String',cbarProps.ylabelstring);
set(get(hndl,'Zlabel'),'String',cbarProps.zlabelstring);

if ~isempty(cbarProps.titleParams)
    set(get(hndl,'Title'),cbarProps.titleParams{:});
    set(hndl,cbarProps.titleParams{:});
end

end

%--------------------------------------------------------------------------

function deleteColorbar()
% Delete the colorbar.
    
% Find the colorbar image.
hImage = findobj(gcf,'tag','CONTOURCMAP','type','image');

if ~isempty(hImage) && ishghandle(hImage)
    % Delete the colorbar axes and image.
    delete(get(hImage(1),'Parent'));
end

end

%-------------------------------------------------------------------------

function setAxesDeleteFcn(ax, axesinfo)

set(ax, 'DeleteFcn', @deleteax)

    function deleteax(~, ~)
        % DeleteFcn callback for the axes holding the colorbar.
        
        if ~isempty(axesinfo) && isfield(axesinfo,'h')
            hDataAxes = axesinfo.h;
            if ishghandle(hDataAxes,'axes');
                % Restore the data axes to its original position.
                set(hDataAxes,'Units','normalized')
                set(hDataAxes,'Position',axesinfo.origpos)
                set(hDataAxes,'Units',axesinfo.units)
            end
        end
    end
end

%-------------------------------------------------------------------------

function deleteim(h, ~)
% DeleteFcn callback for the image object used in the colorbar.

% Double check that h is a 'contourcmap' image, then delete parent axes.
% Let the axes' delete function do the actual work.
if ishghandle(h,'image') && strcmpi(get(h,'Tag'),'contourcmap')
    ax = get(h,'Parent');
    delete(ax)
end

end

%--------------------------------------------------------------------------

function out = canonicalProps(in)
% Expand property names to canonical names

try
    out = lower(validatestring(in, ...
        {'Location', 'ColorAlignment', 'SourceObject', 'TitleString',...
        'XLabelString', 'YLabelString', 'ZLabelString', 'Colorbar'}, ...
        mfilename));
catch e
    if strcmp(e.identifier,'MATLAB:contourcmap:unrecognizedStringChoice') ...
            && ~isempty(e.cause)
        rethrow(e)
    else
        out = in;
    end
end

end
