function h = regionmap(mapFunctionName, varargin)
%REGIONMAP Construct a map axes for a region of the world or USA.

%   Parse and process inputs for either 'worldmap' or 'usamap'.
%   MAPFUNCTIONNAME is the name of one of these functions.  This
%   implementation detects the following obsolete syntaxes, ignores the
%   extraneous input and issuing a warning:
%
%   1) First argument matches a string in {'hi', 'lo', 'allhi'}
%
%   2) Last argument matches a "type" string:
%        {'patch','line','patchonly','lineonly','none','mesh',...
%        'meshonly','mesh3d','dem','demonly','dem3d','dem3donly',..
%        'ldem3d','ldem3donly','lmesh3d','lmesh3donly'}
%
%   3) 'only' is appended to a region or state name

% Copyright 2004-2010 The MathWorks, Inc.
% $Revision: 1.1.6.12.2.1 $  $Date: 2010/12/03 21:43:51 $ 

args = varargin;

% If present, remove obsolete leading argument ('hi', 'lo', or 'allhi')
% and issue warning (worldmap only).
if strcmp(mapFunctionName,'worldmap')
    args = removeObsoleteFirstArg(args);
end

% If present, remove obsolete trailing "type" argument and issue warning.
args = removeObsoleteTypeString(mapFunctionName, args);

% Check argument count now that obsolete args have been removed.
checknargin(0, 2, numel(args), mapFunctionName);

% Set parameters specific to worldmap or usamap
if strcmp(mapFunctionName,'worldmap')
    load('regions','worldRegions');
    validRegions =  worldRegions;
    load('regions','worldAbbreviations')
    abbreviations = worldAbbreviations;
    specialRegions = {};
    constructMapAxes = @constructMapAxesWorld;
    mapLimitIncrement = 5;  % Degrees
else
    load('regions','stateRegions');
    validRegions =  stateRegions;
    load('regions','stateAbbreviations')
    abbreviations = stateAbbreviations;
    specialRegions = {'conus', 'all', 'allequal'};
    constructMapAxes = @constructMapAxesUSA;
    mapLimitIncrement = 1;  % Degrees
end

switch(numel(args))

    case 0    % Syntax: WORLDMAP or USAMAP

        % Ask the user to select a region from a dialog box, returning the
        % result in the string regionName. If the user cancels,
        % regionName will be empty and worldmap (or usamap) is done.
        dialogList = [specialRegions, {validRegions.name}];
        regionName = getRegionFromDialog(dialogList);
        if isempty(regionName)
            h = [];
            return
        else
            % Check for special (USA) map regions.
            if any(strcmpi(regionName,specialRegions))
                % Construct one ('conus') or three ('all','allequal') map
                % axes with hard-coded map projections and limits.
                h = constructSpecialMapAxes(regionName);
            else
                % Lookup the latitude/longitude limits for the region.
                % Embed regionName in a cell array before calling
                % getLimitsFromRegions, because it expects a list of
                % regions.
                [latlim, lonlim] = getLimitsFromRegions(...
                    {regionName}, validRegions, mapLimitIncrement);
                h = constructMapAxes(latlim, lonlim);
            end
        end

    case 1   % Syntax: WORLDMAP REGION or WORLDMAP(REGIONS)
             %         or USAMAP STATE or USAMAP(STATES)

        % If necessary, convert REGIONS from a string or "padded string
        % matrix" to a cell array of strings. If detected, strip off suffix
        % 'only' and issue warning.
        regionNames = checkMapRegions(args{1});

        % Check for special (USA) map regions.
        if isSpecialRegion(mapFunctionName, regionNames, specialRegions)
            % Construct one ('conus') or three ('all','allequal') map axes with
            % hard-coded map projections and limits.
            h = constructSpecialMapAxes(regionNames{1});
        else
            % If necessary, translate from short list of permissible abbreviations.
            if ~isempty(abbreviations)
                regionNames = expandAbbreviations(regionNames, abbreviations);
            end

            % Lookup the latitude/longitude limits for the region and combine
            % limits from multiple regions.
            [latlim, lonlim] = getLimitsFromRegions(...
                regionNames, validRegions, mapLimitIncrement);
            h = constructMapAxes(latlim, lonlim);
        end
    case 2   % Syntax: WORLDMAP(Z, R) or WORLDMAP(LATLIM, LONLIM)
             %         or USAMAP(Z, R) or USAMAP(LATLIM, LONLIM)
             
        s = size(args{2});
        arg2IsRefVectorOrMatrix = isequal(s, [1 3]) || isequal(s, [3 2]);
        if arg2IsRefVectorOrMatrix || isa(args{2}, 'spatialref.GeoRasterReference')
            % WORLDMAP(Z, R)
            [latlim, lonlim] = getLimitsFromDataGrid(args{:}, mapFunctionName);
        else
            % WORLDMAP(LATLIM, LONLIM)
            [latlim, lonlim] = checkMapLimits(args{:}, mapFunctionName);
        end
        h = constructMapAxes(latlim, lonlim);

    otherwise

        error('map:regionmap:internalError', 'More than two input arguments.')
end

% On maps of the whole world, move the meridian labels south a little
% to keep them from overwriting the parallel labels for the equator.
if strcmpi(getm(h(1),'MapProjection'),'robinson')
    setm(h,'MLabelParallel',-10)
end

% Adjust the axes settings to create a more pleasing map figure. 
tightmap

%--------------------------------------------------------------------------

function args = removeObsoleteFirstArg(args)
% If detected, remove an obsolete initial argument from the argument list
% and issue a warning.

obsoleteLeadingArgs = {'hi', 'lo', 'allhi'};

firstArgMatchesObsoleteLeadingArgs = ...
    numel(args) > 0 && ...
    ischar(args{1}) && ...
    any(strcmpi(args{1},obsoleteLeadingArgs));

if firstArgMatchesObsoleteLeadingArgs
    warning('map:regionmap:ignoringFirstArg', ...
        'In function worldmap argument ''%s'' is obsolete and will be ignored.',...
         args{1})
    args(1) = [];
end

%--------------------------------------------------------------------------

function args = removeObsoleteTypeString(mapFunctionName, args)
% If detected, remove an obsolete trailing argument, a "type" string, from
% the argument list and issue a warning.

typeStrings = {'patch','line','patchonly','lineonly','none','mesh',...
               'meshonly','mesh3d','dem','demonly','dem3d','dem3donly',...
               'ldem3d','ldem3donly','lmesh3d','lmesh3donly'};
           
lastArgMatchesTypeString = ...
    numel(args) > 0 && ...
    ischar(args{end}) && ...
    any(strcmpi(args{end},typeStrings));

if lastArgMatchesTypeString
    warning('map:regionmap:ignoringTypeString', ...
        'In function %s argument ''%s'' is obsolete and will be ignored.',...
         mapFunctionName, args{end})
    args(end) = [];
end
   
%--------------------------------------------------------------------------

function region = getRegionFromDialog(regionlist)
% Ask the user to select a region by name from a list dialog.  Return a
% string containing the name of the region.  If the user cancels, return
% empty.

indx = listdlg('ListString',     regionlist,...
               'SelectionMode', 'single',...
               'Name',          'Select a region');

if isempty(indx)
    region = '';
else
    region = regionlist{indx};
end

%--------------------------------------------------------------------------

function result = isSpecialRegion(...
    mapFunctionName, regionNames, specialRegions)
% Check to see if there single region name that matches one of the special
% regions.  However, if multiple regions names are given, then _none_ of
% them may match a special region.

if numel(regionNames) == 1
    result = any(strcmpi(regionNames{1},specialRegions));
elseif numel(regionNames) > 1
    for k = 1:numel(regionNames)
        if any(strcmpi(regionNames{k},specialRegions))
            error('map:regionmap:specialRegionWithOthers', ...
                'The special region, ''%s'', may not be combined with\nother regions or states in function %s.',...
                regionNames{k}, mapFunctionName)
        end
    end
    result = false;
else
    result = false;
end

%--------------------------------------------------------------------------

function regions = expandAbbreviations(regions, abbreviations)
% For each string in REGIONS, translate from a short list of acceptable
% abbreviations if a match is found.  Otherwise, leave the string as-is.

short = abbreviations(:,1);
full  = abbreviations(:,2);
for k = 1:numel(regions)
    index = strcmpi(regions{k}, short);
    if any(index)
        index = find(index);
        regions{k} = full{index(1)};
    end
end

%--------------------------------------------------------------------------

function regionNames = checkMapRegions(regionNames)
% Validate and convert REGIONNAMES to a cell array of strings, removing
% padded blanks and, if present, the obsolete suffix: 'only'. REGIONNAMES
% may be input as a string, a "padded string matrix" (an old syntax that is
% no longer publicized but still supported), or a cell array of strings.

% If regionNames is a cell array, make sure it contains only strings.
if iscell(regionNames)
    for k = 1:numel(regionNames)
        if ~isString(regionNames{k})
            error('map:regionmap:nonStringsInRegion', ...
                'Regions cell array should contain only strings.')
        end
    end
end
  
% If regionNames is a string or "padded string matrix", convert it to a cell array.
if ischar(regionNames)
    regionNames = cellstr(regionNames);
end

% Remove any 'only' suffixes plus trailing or leading spaces.  An 'only'
% suffix is allowed for backward compatibility, but no longer has any
% effect except to trigger a warning.
for k = 1:numel(regionNames)
    regionName = deblank(regionNames{k});
    endsInOnly = (numel(regionName) > 4) && strcmpi(regionName((end-3):end),'only');
    if endsInOnly
        regionName((end-3):end) = [];
        warning('map:regionmap:ignoringSuffixOnly', ...
            'Ignoring suffix ''only'' in region name ''%s''. This feature is obsolete.', ...
            regionNames{k})
    end
    
    % Remove leading blanks and possible trailing blanks
    % (in case suffix was 'only')
    regionNames{k} = strtrim(regionName);
end

%--------------------------------------------------------------------------

function q = isString(s)
% Return true iff input S is a 1-by-N array of class char.

q = ischar(s) && (size(s,1) == 1) && (ndims(s) == 2);

%--------------------------------------------------------------------------

function [latlim, lonlim] = getLimitsFromDataGrid(Z, R, mapFunctionName)
% Compute latitude and longitude limits from a regular data grid and its
% referencing vector or matrix, R.

R = internal.map.convertToGeoRasterRef( ...
    R, size(Z), 'degrees', mapFunctionName, 'R', 2);

latlim = R.Latlim;
lonlim = R.Lonlim;

% Avoid trims if displaying subset of world.
if abs(diff(lonlim)) ~= 360
    lonlim = lonlim + [-1 0]*(10*epsm('deg')); % avoids trimming
end

%--------------------------------------------------------------------------

function [latlim, lonlim] = checkMapLimits(latlim, lonlim, mapFunctionName)
% Check type, size, and values of latitude and longitude limits.

checkinput(latlim, {'double'}, {'real','finite'}, mapFunctionName, 'LATLIM', 1)
checkinput(lonlim, {'double'}, {'real','finite'}, mapFunctionName, 'LONLIM', 2)

if numel(latlim) ~= 2
    error('map:regionmap:latlimWrongLength',...
        'In function %s, %s must have exactly two elements.',...
         mapFunctionName, 'LATLIM')
end

if numel(lonlim) ~= 2
    error('map:regionmap:lonlimWrongLength', ...
         'In function %s, %s must have exactly two elements.',...
          mapFunctionName, 'LONLIM')
end

minlat = -90;  % In degrees
maxlat =  90;  % In degrees
if (latlim(1) < minlat) || (latlim(2) > maxlat)
    error('map:regionmap:latlimOutOfRange', ...
         'In function %s, %s must be within the range [-90 90].',...
          mapFunctionName, 'LATLIM')
end

if latlim(1) >= latlim(2)
    error('map:regionmap:latlimOutOfOrder', ...
         'In function %s, %s must be less than %s.',...
          mapFunctionName, 'LATLIM(1)', 'LATLIM(2)')
end

if lonlim(2) > lonlim(1) + 360
    lonlim(2) = lonlim(1) + 360;
end

%--------------------------------------------------------------------------

function [latlim, lonlim] = getLimitsFromRegions(...
    regionNames, validRegions, inc)
% Return the most compact possible latitude and longitude limits that
% encompass the regions listed by name in cell array regionNames.

latlim = [];
lonlim = [];

validRegionNames = {validRegions.name};

for k = 1:numel(regionNames)
    index = strncmpi(regionNames{k},validRegionNames,numel(regionNames{k}));
    if ~any(index)
        error('map:regionmap:unknownRegion', ...
            'Unknown region name: %s.', regionNames{k})
    end
    index = find(index);
    index = index(1);   % Should be a scalar, but just in case...
    latlim = mergelatlimits(latlim, validRegions(index).latlim);
    lonlim = mergelonlimits(lonlim, validRegions(index).lonlim);
    % Possible enhancements:  Merge all the regions simultaneously to
    % ensure that the longitude limit result is independent of the
    % processing order.  Also ensure that the longitude limits do not span
    % more than 360 degrees.
end

% Snap map limits to increments of INC, with a 1 degree buffer, except
% for limits that are already exact multiples of INC.
buffer = 1;

if mod(latlim(1),inc) ~= 0
    latlim(1) = inc * floor((latlim(1) - buffer)/inc);
end
if mod(lonlim(1),inc) ~= 0
    lonlim(1) = inc * floor((lonlim(1) - buffer)/inc);
end
if mod(latlim(2),inc) ~= 0
    latlim(2) = inc * ceil((latlim(2) + buffer)/inc);
end
if mod(lonlim(2),inc) ~= 0
    lonlim(2) = inc * ceil((lonlim(2) + buffer)/inc);
end

% Ensure that latitude limits remain within [-90 90].
if latlim(1) < -90
    latlim(1) = -90;
end
if latlim(2) > 90
    latlim(2) = 90;
end

%--------------------------------------------------------------------------

function latlim = mergelatlimits(latlim1, latlim2)

% Compute the tightest possible latitude limits encompassing both the
% interval defined by 1-by-2 vector LATLIM1 and the interval defined by
% 1-by-2 vector in LATLIM2.  Note that either input could be empty.
limits = [latlim1 latlim2];
latlim = [min(limits) max(limits)];

%--------------------------------------------------------------------------

function lonlim = mergelonlimits(lonlim1, lonlim2)

% Compute the tightest possible longitude limits encompassing both the
% interval defined by 1-by-2 vector LONLIM1 and the interval defined by
% 1-by-2 vector LONLIM2.  In addition, LONLIM1, LONLIM2, or both may be
% empty.

if isempty(lonlim1)
    lonlim = lonlim2;
elseif isempty(lonlim2)
    lonlim = lonlim1;
else
    % Shift both intervals such that the first one starts at zero.  Call
    % the shifted versions i1 and i2.
    s1 = lonlim1(1);
    i1 = zero22pi(lonlim1 - s1);
    i2 = lonlim2 - s1;
    if zero22pi(i2(1)) <= i1(2)
        % We have overlap, with interval 2 starting within interval 1
        s2 = i2(1) - zero22pi(i2(1));
        % If necessary, shift i2 by a multiple of 360 degrees, so that
        % i2(1) falls within i1.  Call the result j2.
        j2 = i2 - s2;
        % Merge i1 and j2
        j = [0, max(i1(2),j2(2))];
    elseif zero22pi(i2(2)) <= i1(2)
        % We have overlap, with interval 2 ending within interval 1
        s2 = i2(2) - zero22pi(i2(2));
        % If necessary, shift i2 by a multiple of 360 degrees, so that
        % i2(2) falls within i1.  Call the result j2.
        j2 = i2 - s2;
        % Merge i1 and j2
        j = [min(0,j2(1)) i1(2)];
    else
        % Neither overlap condition was met; there is no overlap. We can
        % define j (shifted output interval) by either putting i2 to the
        % east of i1, or by putting it to the west.  We'll make the choice
        % that minimizes the width of j.
        width1 = zero22pi(i2(2));   % Width putting i2 to the east
        width2 = i1(2) - (zero22pi(i2(1)) - 360);  % Width putting i2 to the west
        if width1 <= width2
            j = [0 width1];
        else
            j = i1(2) + [-width2 0];
        end
    end
    % Undo the shift s1
    lonlim = j + s1;
end

% The following changes the behavior of usamap('alaska') and seems
% unnecessary:
%
% lonlim = npi2pi(lonlim);  Return results in the interval [-180 180].

%--------------------------------------------------------------------------

function ax = constructMapAxesWorld(latlim,lonlim)

% Construct a map axes suitable to the specified latitude and longitude limits.

% Ensure row vectors.
latlim = latlim(:)';
lonlim = lonlim(:)';

% Ensure ascending lonlim.
if lonlim(1) > lonlim(2)
    if lonlim(1) > 180
        lonlim(1) = lonlim(1) - 360;
    else
        lonlim(2) = lonlim(2) + 360;
    end
end

% Compute a nice increment for labels and grid.
% Pick a ticksize that gives at least 3 grid lines

mindiff = min(abs([diff(latlim) diff(lonlim)]));
ticcandidates = [.1 .5 1 2 2.5 5 10 15 20 30 45 ] ;
[~,indx] = min( abs( 3 - (mindiff ./ ticcandidates) ));

ticksize = ticcandidates(indx);
roundat = 0;
if mod(ticksize,1) ~= 0; roundat = -1; end

% Select a projection based on latlim,lonlim.

if isequal(latlim,[-90 90]) && diff(lonlim) >= 360 % entire globe
   projection = 'robinson';
elseif max(abs(latlim)) < 20% straddles equator, but doesn't extend into extreme latitudes
   projection = 'mercator';
elseif abs(diff(latlim)) <= 90 && abs(sum(latlim)) > 20 && max(abs(latlim)) < 90 % doesn't extend to the pole, not stradling equator
   projection = 'eqdconic';
elseif abs(diff(latlim)) < 85 && max(abs(latlim)) < 90 % doesn't extend to the pole, not stradling equator
   projection = 'sinusoid';
elseif max(latlim) == 90 && min(latlim) >= 84
   projection = 'upsnorth';
elseif min(latlim) == -90 && max(latlim) <= -80
   projection = 'upssouth';
elseif max(abs(latlim)) == 90 && abs(diff(lonlim)) < 180
   projection = 'polycon';
elseif max(abs(latlim)) == 90 
   projection = 'eqdazim';
else
   projection = 'miller';
end

% Delete existing map axes if necessary
if ismap(gca) == 1
   delete(get(gca,'Children'));
end

% Set up the axes with the selected projection and parameters.
% More than one path because of error messages setting parallels
% when projections don't have any.

if strcmp(projection,'upsnorth') || strcmp(projection,'upssouth')
   
   ax = axesm('MapProjection','ups','Zone',projection(4:end),...
      'Frame','on',...
      'Grid','on',...
      'LabelRotation','on',...
      'MeridianLabel','on',...
      'ParallelLabel','on',...
      'MLabelParallel',0 ...
      );
   
else
   
   mstruct = defaultm(projection);
   
   if strcmp(projection,'eqdazim')
       % Separate graticule meridians and meridian labels by 30 degrees for
       % polar azimuthal projections.
       lonGratSpacingAzimuthal = 30;
       
       ax = axesm(...
         'MapProjection',projection,...
         'FLatLimit',[-Inf abs(diff(latlim))],...
         'Origin',[90*sign(latlim(2)) mean(lonlim) 0],...
         'MLabelLoc',lonGratSpacingAzimuthal,...
         'MLineLoc',lonGratSpacingAzimuthal,...
         'PLabelLoc',ticksize,...
         'PLineLoc',ticksize,...
         'MLabelRound',roundat,...
         'PLabelRound',roundat,...
         'GColor',[.75 .75 .75],...
         'GLineStyle',':',...
         'Frame','on',...
         'Grid','on',...
         'LabelRotation','on',...
         'MeridianLabel','on',...
         'ParallelLabel','on',...
         'MLabelParallel',0 ...
         );
      
   elseif mstruct.nparallels > 0
         
       ax = axesm(...
         'MapProjection',projection,...
         'MapLatLimit',latlim,...
         'MapLonLimit',lonlim,...
         'MapParallels',[],...
         'MLabelLoc',ticksize,...
         'MLineLoc',ticksize,...
         'PLabelLoc',ticksize,...
         'PLineLoc',ticksize,...
         'MLabelRound',roundat,...
         'PLabelRound',roundat,...
         'MLabelPar',0,...
         'GColor',[.75 .75 .75],...
         'GLineStyle',':',...
         'Frame','on',...
         'Grid','on',...
         'LabelRotation','on',...
         'MeridianLabel','on',...
         'ParallelLabel','on'...
         );
   else
      ax = axesm(...
         'MapProjection',projection,...
         'MapLatLimit',latlim,...
         'MapLonLimit',lonlim,...
         'MLabelLoc',ticksize,...
         'MLineLoc',ticksize,...
         'PLabelLoc',ticksize,...
         'PLineLoc',ticksize,...
         'MLabelRound',roundat,...
         'PLabelRound',roundat,...
         'MLabelPar',0,...
         'GColor',[.75 .75 .75],... 
         'GLineStyle',':',...
         'Frame','on',...
         'Grid','on',...
         'LabelRotation','on',...
         'MeridianLabel','on',...
         'ParallelLabel','on'...
         );
   end
end

set(ax, 'Visible', 'off')
set(get(ax,'Title'), 'Visible', 'on')

%--------------------------------------------------------------------------

function h = constructMapAxesUSA(latlim,lonlim)
% Construct a map axes suitable to the specified latitude and longitude limits,
% for maps covering all or part of the Conterminous U.S.

% Ensure row vectors.
latlim = latlim(:)';
lonlim = lonlim(:)';

% Ensure ascending lonlim.
if lonlim(1) > lonlim(2)
    if lonlim(1) > 180
        lonlim(1) = lonlim(1) - 360;
    else
        lonlim(2) = lonlim(2) + 360;
    end
end

% Compute a nice increment for labels and grid.
% Pick a ticksize that gives at least 3 grid lines

mindiff = min(abs([diff(latlim) diff(lonlim)]));
ticcandidates = [.1 .5 1 2 2.5 5 10 15 20 30 45 ] ;
[~,indx] = min( abs( 4 - (mindiff ./ ticcandidates) ));

ticksize = ticcandidates(indx);
roundat = 0;
if mod(ticksize,1) ~= 0; roundat = -1; end

% States are small enough to display conformally with little error
projection = 'lambert';

% set up the map axes
h = axesm(...
    'MapProjection',projection,...
    'MapLatLimit',latlim,...
    'MapLonLimit',lonlim,...
    'MapParallels',[],...
    'MLabelLoc',ticksize,...
    'MLineLoc',ticksize,...
    'PLabelLoc',ticksize,...
    'PLineLoc',ticksize,...
    'MLabelRound',roundat,...
    'PLabelRound',roundat,...
    'MLabelPar',0,...
    'GColor',[.75 .75 .75],...
    'GLineStyle',':',...
    'Frame','on',...
    'Grid','on',...
    'LabelRotation','on',...
    'MeridianLabel','on',...
    'ParallelLabel','on'...
    );

set(h, 'Visible', 'off')
set(get(h,'Title'), 'Visible', 'on')

% Leave a little space around non-cylindrical projections to avoid clipping the frame
names = maps('idlist');
pclass = maps('classcodes');
indx = strcmp(projection,names);
thispclass = deblank(pclass(indx,:));

switch thispclass
case 'Cyln'
   tightmap tight
otherwise
   tightmap loose
end

%--------------------------------------------------------------------------

function h = constructSpecialMapAxes(region)
% Construct map axes for the United States for three special
% region designations:  'all', 'allequal', or 'conus'.

% Construct either one or three axes, with carefully-defined positions.
switch(lower(region))
    case 'conus'
        conusPosition = [ -0.0078125    0.13542     1.0156    0.73177 ];
        h = constructConusMapAxes(conusPosition);
        
    case 'all'
        conusPosition  = [ -0.0078125    0.13542     1.0156    0.73177 ];
        alaskaPosition = [  0.029297     0.10938     0.39453   0.36458 ];
        hawaiiPosition = [  0.31055      0.17188     0.19141   0.20052 ];
        hConus =  constructConusMapAxes(conusPosition);
        hAlaska = constructAlaskaMapAxes(alaskaPosition);
        hHawaii = constructHawaiiMapAxes(hawaiiPosition);
        h = [hConus hAlaska hHawaii];
        set(gcf,'CurrentAxes',hConus)
        
    case 'allequal'
        conusPosition  = [0.1        0.25     0.85 0.6];
        alaskaPosition = [0.019531  -0.020833  0.2 0.2];
        hawaiiPosition = [0.6        0.2       0.2 0.2];
        hConus =  constructConusMapAxes(conusPosition);
        hAlaska = constructAlaskaMapAxes(alaskaPosition);
        hHawaii = constructHawaiiMapAxes(hawaiiPosition);
        h = [hConus hAlaska hHawaii];
        axesscale(hConus) % Resize axes for Alaska and Hawaii.
        set(gcf,'CurrentAxes',hConus)
end

% Hide the axes, but not the title of the first (CONUS) axes.
set(h, 'Visible', 'off')
set(get(h(1),'Title'), 'Visible', 'on')

%--------------------------------------------------------------------------

function hConus = constructConusMapAxes(position)
% Map axes for CONUS.

hConus = axes('Position', position);

axesm(...
    'MapProjection','eqaconic',...
    'MapParallels', [29.5 49.5],...
    'Origin', [0 -100 0],...
    'FLatLimit', [23 50],...
    'FLonLimit', [-27  35],...
    'Frame','on',...
    'Grid','on',...
    'MeridianLabel','on',...
    'ParallelLabel','on',...
    'MLineLocation', 10, 'PLineLocation', 10,...
    'MLabelLocation', 10, 'PLabelLocation', 10,...
    'LabelRotation', 'on',...
    'MLabelPar', 0)

xlim([-0.56 0.6])   % Put nice amount of white space around CONUS, and lock down
ylim([ 0.30 0.92])	% limits to disable autoscaling. Use AXIS AUTO to undo.

%--------------------------------------------------------------------------

function hAlaska = constructAlaskaMapAxes(position)
% Map axes for Alaska. Scale is considerably different from CONUS.
        
hAlaska = axes('Position', position);

axesm(...
    'MapProjection','eqaconic',...
    'MapParallels', [55 65],...
    'Origin', [0 -150 0],...
    'FLatLimit', [50 75],...
    'FLonLimit', [-25  25],...
    'Frame','on',...
    'Grid','on',...
    'MeridianLabel','on',...
    'ParallelLabel','on',...
    'MLineLocation', 10, 'PLineLocation', 10,...
    'MLabelLocation', 10, 'PLabelLocation', 10,...
    'LabelRotation', 'on',...
    'MLabelPar', 0)

xlim([-0.3 0.3 ])
ylim([ 0.7  1.3])

%--------------------------------------------------------------------------

function hHawaii = constructHawaiiMapAxes(position)
% Map axes for Hawaii.

hHawaii = axes('Position',position);

axesm(...
    'MapProjection','eqaconic',...
    'MapParallels', [8 18],...
    'Origin', [0 -157 0],...
    'FLatLimit', [18 24],...
    'FLonLimit', [-3  3],...
    'Frame','on',...
    'Grid','on',...
    'MeridianLabel','on',...
    'ParallelLabel','on',...
    'LabelRotation', 'on',...
    'MLineLocation', 1, 'PLineLocation', 1,...
    'MLabelLocation', 5, 'PLabelLocation', 5)

xlim([-0.064103  0.058974])
ylim([ 0.305190  0.440260])
