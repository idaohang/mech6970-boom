function [Z, R, levelList, pvpairs] = parseContourInputs(args)
%PARSECONTOURINPUTS Parse inputs for contourm
%
%   [Z, R, levelList, pvpairs] = parseContourInputs(ARGS) parses the
%   inputs to the CONTOURM function.  The output PVPAIRS is a property
%   name-value cell array in which the names have been expanded and
%   verified to be valid CONTOURM property names.

% Copyright 2004-2010 The MathWorks, Inc.
% $Revision: 1.1.6.12 $  $Date: 2010/11/17 11:25:14 $

% Get data grid and referencing object from arguments 1,2 or 1,2,3.
if (numel(args) == 2) ...
        || isequal(size(args{2}),[1 3]) ...
        || isequal(size(args{2}),[3 2]) ...
        || isobject(args{2})
    % Using (Z, R, ...) syntax.
    
    Z = args{1};
    R = args{2};
    
    args(1:2) = [];
    
    assert(~isempty(Z), 'map:validate:expectedNonEmpty', ...
        'Input data %s must not be empty.', 'Z')
    
    %  If R is already spatial referencing object, validate it. Otherwise
    %  convert the input referencing vector or matrix.
    R = internal.map.convertToGeoRasterRef( ...
        R, size(Z), 'degrees', 'contourm', 'R', 2);
else
    % Using (LAT, LON, Z, ...) syntax.
    
    lat = args{1};
    lon = args{2};
    Z   = args{3};
    
    args(1:3) = [];
    
    assert(~isempty(Z), 'map:validate:expectedNonEmpty', ...
        'Input data %s must not be empty.', 'Z')

    % Validate sizes and consistency of latitude, longitude, and data
    % arrays; expand latlim/lonlim inputs.
    R = checkLatLonZ(lat, lon, Z);    
end

% Warn if Z is uniform or contains no finite values.
checkZ(Z)

% Extract the levels argument (N or V) if provided.
[levelsArg, args] = checkForLevelsArg(args);

% Decode the linespec argument, if provided.
[lineprops, args] = checkLinespec(args);

% Expand and validate property in the list of parameter-value pairs, if
% provided. (And after calling checkLinespec, numel(args) is even.)
pvpairs = checkParameterNames([lineprops args]);

% Get the level list vector if provided, construct if not.
% If pvpairs contains a LevelStep pair, remove it.
[levelList, pvpairs] = checkLevels(Z, levelsArg, pvpairs);

%-----------------------------------------------------------------------

function checkZ(Z)
% Check for non-finite or uniform input grid Z.

k = find(isfinite(Z));
if ~any(k)
    warning('map:contourm:nonFiniteData', ...
        'Contours not rendered for non-finite %s.','Z');
else
    if min(Z(k)) == max(Z(k))
        warning('map:contourm:uniformData', ...
            'Contour lines not rendered for uniform %s.','Z');
    end
end

%-----------------------------------------------------------------------

function [lineprops, args] = checkLinespec(args)

lineprops = {};
if mod(numel(args),2) == 1
    linespec = args{1};
    args(1) = [];
    [lineStyle, lineColor, marker, msg] = colstyle(linespec);
    if ~isempty(msg)
        error('map:contourm:invalidColorSpec', ...
            'The input ''%s'' is not a valid color specification.', ...
            linespec)
    end
    if ~isempty(marker)
        warning('map:contourm:ignoringMarker', ...
            'Contour lines do not have markers; marker %s will be ignored.', ...
            marker)
    end
    if ~isempty(lineStyle)
        lineprops = {'LineStyle',lineStyle};
    end
    if ~isempty(lineColor)
        lineprops(end+1:end+2) = {'LineColor',lineColor};
    end
end

%-----------------------------------------------------------------------

function R = checkLatLonZ(lat, lon, Z)

validateattributes(lat, {'double'}, {'2d','real','finite'}, 'contourm', 'LAT')
validateattributes(lon, {'double'}, {'2d','real','finite'}, 'contourm', 'LON')
validateattributes(lat, {'double'}, {'2d','real','finite'}, 'contourm', 'Z')

if numel(lat) == 2 && numel(lon) == 2 && ~isequal(size(Z), size(lat))
    % latlim, lonlim syntax
    R = georasterref('RasterSize',size(Z), ...
        'Latlim', lat, 'Lonlim', lon, 'RasterInterpretation', 'cells');
else
    % Check Latitude.
    if ~isequal(size(lat), size(Z))
        if min(size(lat))==1 && length(lat) ~= size(Z,1)
            error('map:contourm:invalidLatLength', ...
                'Length of %s must match number of rows in %s.', 'LAT', 'Z')
        elseif min(size(lat)) ~= 1
            error('map:validate:inconsistentSizes', ...
                '%s and %s must have the same size.', 'LAT', 'Z')
        end
    end
    
    % Check Longitude.
    if ~isequal(size(lon), size(Z))
        if min(size(lon))==1 && length(lon) ~= size(Z,2)
            error('map:contourm:invalidLonLength', ...
                'Length of %s must match number of columns in %s.', 'LON', 'Z')
        elseif min(size(lon))~=1
            error('map:validate:inconsistentSizes', ...
                '%s and %s must have the same size.', 'LON', 'Z')
        end
    end
    
    if isvector(lat) && isvector(lon)
        % Construct mesh if required.
        [lon, lat] = meshgrid(lon, lat);
    end
    
    % Construct a "geomesh" structure.
    R = struct('AngleUnits','degrees','LatMesh',lat,'LonMesh',lon);
end

%-----------------------------------------------------------------------

function [levelsArg, args] = checkForLevelsArg(args)

if ~isempty(args) && ~ischar(args{1})
    % (LAT, LON, Z, N, ...) or (Z, R, N) or
    % (LAT, LON, Z, LEVELS, ...) or (Z, R, LEVELS)
    levelsArg = args{1};
    args(1) = [];
    validateattributes(levelsArg, {'numeric'},{'real','finite','vector'}, ...
        'contourm','N or V')
else
    % N and LEVELS were both omitted.
    levelsArg = [];
end

%-----------------------------------------------------------------------

function [levelList, pvpairs] = checkLevels(Z, levelsArg, pvpairs)
% Input levelsArg can be:
%
%   N  -- The number of contour levels to use
%   V  -- A levelList vector
%   [] -- If the arguments N and LEVELS were omitted

Z = real(Z);
k = find(isfinite(Z));
zmax = max(Z(k));
zmin = min(Z(k));

[levelStep, pvpairs] = extractLevelStepParameter(pvpairs);
[levelList, pvpairs] = extractLevelListParameter(pvpairs);

if isempty(levelList)
    % Consider N, V, and LevelStep only if LevelList was not specified.
    if ~isempty(levelsArg)
        % (LAT, LON, Z, N, ...) or (Z, R, N) or
        % (LAT, LON, Z, LEVELS, ...) or (Z, R, LEVELS)
        
        % Adapted from: toolbox/matlab/specgraph/contour.m (parseargs)
        if isscalar(levelsArg) && (levelsArg > 0) && (fix(levelsArg) == levelsArg)
            % N -- A positive integer value
            if isempty(zmax)
                % No finite data
                levelList = 0;
            elseif levelsArg == 1
                levelList = (zmin + zmax)/2;
            else
                levelList = linspace(zmin, zmax, levelsArg + 2);
                levelList = levelList(2:end-1);
            end
            
            if ~isempty(levelStep)
                % Let user know that LevelStep parameter is being ignored.
                warning('map:contourm:ignoringLevelStep1', ...
                    ['Ignoring %s parameter because the number of', ...
                    ' contour levels was also specified.'], 'LevelStep')
            end
        else
            % LEVELS
            levelList = unique(levelsArg);
            
            if ~isempty(levelStep)
                % Let user know that LevelStep parameter is being ignored.
                warning('map:contourm:ignoringLevelStep2', ...
                    ['Ignoring %s parameter because a contour levels', ...
                    ' vector was also specified.'], 'LevelStep')
            end
        end
    else
        if isempty(levelStep)
            % Neither the number of levels nor a level list are specified,
            % and the LevelStep parameter was not supplied:
            % select a level step using a heuristic based on the data range.
            levelStep = selectLevelStep(zmin, zmax);
        end
        
        levelList = constructLevelList(zmin, zmax, levelStep);
    end
else
    if ~isempty(levelStep)
        % Let user know that LevelStep parameter is being ignored.
        warning('map:contourm:ignoringLevelStep3', ...
            ['Ignoring %s parameter because a list of', ...
            ' contour levels was also specified.'], 'LevelStep')
    end
end

%-----------------------------------------------------------------------

function [levelstep, pvpairs]  = extractLevelStepParameter(pvpairs)
% Look (carefully) for a level step value in the pvpairs array.

k = find(strcmp('LevelStep', pvpairs));
if ~isempty(k) && (k(1) < numel(pvpairs)) && (mod(k(1),2) == 1)
    k = k(1);
    levelstep = pvpairs{k+1};
    pvpairs([k k+1]) = [];
    try
    validateattributes(levelstep, {'double'}, ...
        {'positive','real','finite','scalar'}, ...
        'contourm', 'the value of LevelStep')
    catch exception
        warning('map:contourm:ignoringInvalidLevelStep', ...
            'Ignoring %s: %s', 'LevelStep', exception.message')
        levelstep = [];
    end
else
    levelstep = [];
end

%-----------------------------------------------------------------------

function [levelList, pvpairs]  = extractLevelListParameter(pvpairs)
% Look (carefully) for a level list value in the pvpairs array.

k = find(strcmp('LevelList', pvpairs));
if ~isempty(k) && (k(1) < numel(pvpairs)) && (mod(k(1),2) == 1)
    k = k(1);
    levelList = pvpairs{k+1};
    pvpairs([k k+1]) = [];
    try
    validateattributes(levelList, {'double'}, ...
        {'real','finite','vector'}, ...
        'contourm', 'the value of LevelList')
    catch exception
        warning('map:contourm:ignoringInvalidLevelList', ...
            'Ignoring %s: %s', 'LevelList', exception.message')
        levelList = [];
    end
else
    levelList = [];
end

%-----------------------------------------------------------------------

function levelstep = selectLevelStep(zmin, zmax)
% Select a level step automatically.

% Adapted from:
%    toolbox/matlab/specgraph/@specgraph/@contourgroup/refresh.m

range = zmax-zmin;
range10 = 10^(floor(log10(range)));
nsteps = range/range10;
if nsteps < 1.2
    range10 = range10/10;
elseif nsteps < 2.4
    range10 = range10/5;
elseif nsteps < 6
    range10 = range10/2;
end
levelstep = range10;

%-----------------------------------------------------------------------

function levelList = constructLevelList(zmin, zmax, levelstep)

% Adapted from:
%    toolbox/matlab/specgraph/@specgraph/@contourgroup/refresh.m

if zmin < 0 && zmax > 0
    step = levelstep;
    neg = -step:-step:zmin;
    pos = 0:step:zmax;
    levelList = [fliplr(neg) pos];
elseif zmin < 0 && zmin < zmax
    step = levelstep;
    start = zmin - (step - mod(-zmin,step));
    levelList = start+step:step:zmax;
elseif 0 <= zmin && zmin < zmax
    step = levelstep;
    start = zmin + (step - mod(zmin,step));
    levelList = start:step:zmax;
else
    % zmin == zmax
    levelList = zmin;
end

%-----------------------------------------------------------------------

function args = checkParameterNames(args)
% Set and check the parameter name-value pairs. Assume that args has an
% even number of elements. Allow the names 'EdgeColor' and 'Color', but
% convert them to 'LineColor'. Allow 'DefaultLineColor' and 'LineZ' for
% internal use.

% Names that are valid for users of contourm, contour3m, contourfm.
validNames = sort({ ...
    'LevelStep', ...    % Special parameter; converted to level list
    'LevelList', ...    % Special parameter (level list specified directly)
    'ShowText', ...     % Converted to ContourLines property
    'Fill', ...         % GeoContourGroup property
    'LabelSpacing', ... % GeoContourGroup property
    'LineColor', ...    % GeoContourGroup property
    'LineStyle', ...    % GeoContourGroup property
    'LineWidth', ...    % GeoContourGroup property
    'HandleVisibility', ... % hggroup property
    'Parent', ...           % hggroup property
    'Tag', ...              % hggroup property
    'UserData', ...         % hggroup property
    'Visible'});            % hggroup property

for k = 1:2:numel(args)
    try
        % Full set of allowable names.
        name = validatestring(args{k}, [validNames ...
            {'Color','EdgeColor','DefaultLineColor','LineZ','FillZ'}], ...
            'contourm');
    catch exception
        if strcmp(exception.identifier, ...
                'MATLAB:contourm:unrecognizedStringChoice')
            % Reproduce the error, but include only user-visible names
            % in the message, and change the function name reported in
            % the message from 'contourm' to 'contour3m' or 'contourfm'
            % if either of those fall immediately above 'contourm' in
            % the stack. (Assume that parseContourInputs was called by
            % contourm.)
            functionName = 'contourm';
            namesOnStack = {exception.stack.name};
            if ismember('contour3m',namesOnStack)
                functionName = 'contour3m';
            elseif ismember('contourfm',namesOnStack)
                functionName = 'contourfm';
            end
            validatestring(args{k}, validNames, functionName)
        else
            % Something unexpected happened.
            rethrow(e)
        end
    end
    
    if strcmp(name,'EdgeColor')
        name = 'LineColor';
    elseif strcmp(name,'Color')
        name = 'LineColor';
    end
    args{k} = name;
end

% If 'LineColor' is unspecified, but 'DefaultLineColor' is provided,
% change 'DefaultLineColor' to 'LineColor'. If both are provided, remove
% the 'DefaultLineColor' pair.
if ~any(strcmp('LineColor',args(1:2:end)))
    k = find(strcmp('DefaultLineColor',args(1:2:end)));
    args(1 + 2*(k-1)) = {'LineColor'};
else
    k = find(strcmp('DefaultLineColor',args(1:2:end)));
    args([1 + 2*(k-1), 2*k]) = [];
end

% Filter the LineColor value of 'auto', converting it to 'flat'.
for k = find(strcmp('LineColor',args(1:2:end)));
    j = 2*k;
    value = args{j};
    if strncmpi(value, 'auto', numel(value))
        args{j} = 'flat';
    end
end

% Filter the ShowText parameter, converting to the ContourLabels
% property of internal.mapgraph.GeoContourGroup.
for k = find(strcmp('ShowText',args(1:2:end)));
    j = 2*k;
    args{j-1} = 'ContourLabels';
    value = args{j};
    value = validatestring(value, {'off','on'}, 'contourm', ...
        'value of ShowText property');
    if strcmp(value,'on')
        args{j} = 'all';
    else
        args{j} = 'none';
    end
end
