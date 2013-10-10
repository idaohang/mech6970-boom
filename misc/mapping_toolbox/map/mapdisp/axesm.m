function [h,msg] = axesm(varargin)
%AXESM Define map axes and set map properties
%
%  AXESM activates a GUI to define a map projection for the current axes.
%
%  AXESM(PROPERTYNAME, PROPERTYVALUE,...) uses the map properties in the
%  input list to define a map projection for the current axes. For a list
%  of map projection properties, execute GETM AXES.  All standard
%  (non-mapping) axes properties are controlled using the axes command. For
%  a list of available projections, execute MAPS.
%
%  AXESM(PROJFCN,...) uses the function specified by the string PROJFCN to
%  initialize the map projection.  If this form is used, PROJFCN must
%  correspond to a function on the user's path and is case-sensitive.
%
%  See also AXES, GETM, MAPLIST, MAPS, MFWDTRAN, MINVTRAN, PROJFWD,
%  PROJINV, PROJLIST, SETM.

% Copyright 1996-2010 The MathWorks, Inc.
% $Revision: 1.14.4.19 $  $Date: 2010/03/04 16:22:34 $

% The following syntax is invoked from within Mapping Toolbox(TM)
% function SETM, but is not intended for general use.
%
%    AXESM(MSTRUCT,...) uses the map structure specified by MSTRUCT to
%                       initialize the map projection.

% Obsolete syntax
% ---------------
% [h,msg] = AXESM(...) returns a string indicating any error encountered.
if nargout > 1
    warnObsoleteMSGSyntax(mfilename)
    msg = '';
end

%  Initialize output variables
if nargout ~= 0
    h = [];
end

%  Initialize default map structure.
%  Save each of the field names to compare with passed in properties

mstruct = initmstruct;            %  AXESM algorithm requires mfields to
mfields = fieldnames(mstruct);    %  always remain a cell array.

%  Test input arguments

if (nargin > 0) && any(ishghandle(varargin{1}))
    error(['map:' mfilename ':invalidFirstParam'], ...
        'First argument must be a map property name or map structure')
end

if nargin == 0
    % AXESM
    if ~ismap(gca)
        [~,defproj] = maplist; % get the default projection
        mstruct.mapprojection = defproj;
        mstruct = feval(mstruct.mapprojection,mstruct);
        mstruct = resetmstruct(mstruct);
        set(gca,...
            'NextPlot','Add', ...
            'UserData',mstruct, ...
            'DataAspectRatio',[1 1 1], ...
            'Box','on',...
            'ButtonDownFcn','uimaptbx')
        %  May not hit mapprojection case with non-map axes
        set(gca,'XLimMode','auto','YLimMode','auto')
        showaxes off
    end
    cancelflag = axesmui;
    if nargout ~= 0
        h = cancelflag;
    end
elseif nargin == 1 && ~ischar(varargin{1}) && ~isstruct(varargin{1})
    gcm(varargin{1})
    axes(varargin{1}); %#ok<MAXES>
else
    if rem(nargin,2)
        if isstruct(varargin{1})
            % AXESM(MSTRUCT,...)
            mstruct   = varargin{1};
            startpt   = 2;
            newfields  = sortrows(char(fieldnames(mstruct)));
            testfields = sortrows(char(mfields));
            if any(size(testfields) ~= size(newfields)) ...
                    || any(any(testfields ~= newfields))
                error('map:axesm:invalidMap', ...
                    'Incorrect map structure supplied');
            end
        else
            % AXESM(PROJFCN,...)
            startpt = 2;
            try
                mstruct.mapprojection = maps(varargin{1});
                mstruct = feval(mstruct.mapprojection,mstruct);
            catch %#ok<CTCH>
                if exist(varargin{1},'file') == 2
                    mstruct = feval(varargin{1},mstruct);
                else
                    error('map:axesm:undefinedMapAxis', ...
                        'Unspecified map axis property value')
                end
            end
        end
    else
        % AXESM(PROPERTYNAME, PROPERTYVALUE,...)
        startpt = 1;
    end

    %  Permute the property list so that 'angleunits' (if supplied) is first
    %  and 'mapprojection' is second.  'angleunits' must be processed first
    %  so that all defaults end up in the proper units.  'mappprojection'
    %  must be processed second so that any supplied parameter, such as
    %  'parallels', is not overwritten by the defaults in the projection
    %  function.
    varargin(1:(startpt-1)) = [];
    varargin = reorderprops(varargin);

    % Assign variables to test if a reset property is being set to 'reset'.
    resetFrameGrat = false;
    resetProperties = {'grid', 'meridianlabel', 'parallellabel' ,'frame'};
    
    %  Cycle through the property name-value pairs.
    for j = 1 : 2 : numel(varargin)
        [propname, propvalue] = validateprop(varargin, mfields, j);
        if isequal(propvalue, 'reset') && ...
                any(strncmp(propname, resetProperties, numel(propname)))
             resetFrameGrat = true;
        end
        mstruct = setprop(mstruct, propname, propvalue);
    end

    %  Remove possible NaN left by setMapLatLimit.
    mstruct.origin(isnan(mstruct.origin)) = 0;
    
    %  Check for defaults to be computed.
    mstruct = resetmstruct(mstruct);
    
    if resetFrameGrat
        % Update the axes with the new mstruct.
        set(gca, 'UserData', mstruct);
    else
        % Set GCA to be a map axes
        setgca(mstruct)
    end

    %  Display map frame, lat-lon graticule, and meridian and parallel
    %  labels, if necessary
    setframegrat(mstruct)

    %  Set output variable if necessary
    if nargout >= 1
        h = gca;
    end
end

%-----------------------------------------------------------------------

function props = reorderprops(props)
% Permute the property list, moving the following properties (when
% present) to the start and ordering them as listed here: angleunits,
% mapprojection, zone, origin, flatlimit, flonlimit, maplatlimit,
% maplonlimit.  Also, convert all the property names to lower case.

% Reshape to a 2-by-N: Property names are in row 1, values in row 2.
props = reshape(props,[2,numel(props)/2]);

% Convert all property names to lower case.
props(1,:) = lower(props(1,:));

% Index properties: 101, 102, 103, ...
indx = 100 + (1:size(props,2));

% Determine which columns contain 'angleunits', 'mapprojection', etc.
indx(strmatch('an',    props(1,:))) = 1;  % 'angleunits'
indx(strmatch('mappr', props(1,:))) = 2;  % 'mapprojection'
indx(strmatch('z',     props(1,:))) = 3;  % 'zone'
indx(strmatch('o',     props(1,:))) = 4;  % 'origin'
indx(strmatch('fla',   props(1,:))) = 5;  % 'flatlimit'
indx(strmatch('flo',   props(1,:))) = 6;  % 'flonlimit'
indx(strmatch('mapla', props(1,:))) = 7;  % 'maplatlimit'
indx(strmatch('maplo', props(1,:))) = 8;  % 'maplonlimit'

% Sort indx and save the required permutations in indexsort.
[~, indexsort] = sort(indx);

% Permute the columns of props.
props = props(:,indexsort);

% Turn props back into a row vector.
props = props(:)';

%-----------------------------------------------------------------------

function [propname, propvalu] = validateprop(props, mfields, j)

%  Get the property name and test for validity.
indx = strmatch(lower(props{j}),mfields);
if isempty(indx)
    error('map:axesm:unknownProperty', ...
        'Unrecognized property:  %s',props{j})
elseif length(indx) == 1
    propname = mfields{indx};
else
    error('map:axesm:nonUniqueName', ...
        'Property %s name not unique - supply more characters', ...
        props{j})
end

%  Get the property value, ensure that it's a row vector and convert
%  string-valued property values to lower case.
propvalu = props{j+1};
propvalu = propvalu(:)';
if ischar(propvalu)
    propvalu = lower(propvalu);
end

%-----------------------------------------------------------------------

function mstruct = setprop(mstruct, propname, propvalu)

switch propname

    %*************************************
    %  Properties That Get Processed First
    %*************************************

    case 'angleunits'
        mstruct = setAngleUnits(mstruct,propvalu);

    case 'mapprojection'
        mstruct = setMapProjection(mstruct, propvalu);
        
    case 'zone'
        mstruct = setZone(mstruct, propvalu);

    case 'origin'
        mstruct = setOrigin(mstruct, propvalu);

    case 'flatlimit'
        mstruct = setflatlimit(mstruct, propvalu);

    case 'flonlimit'
        mstruct = setflonlimit(mstruct, propvalu);

    case 'maplatlimit'
        if isempty(propvalu)
            mstruct.maplatlimit = [];
        else
            mstruct = setMapLatLimit(mstruct, propvalu);
        end

    case 'maplonlimit'
        if isempty(propvalu)
            mstruct.maplonlimit = [];
        else
            mstruct = setMapLonLimit(mstruct, propvalu);
        end

    %************************
    %  General Map Properties
    %************************

    case 'aspect'
        validparm = {'normal','transverse'};
        indx = strmatch(propvalu,validparm);
        if length(indx) == 1
            mstruct.aspect = validparm{indx};
        else
            error(['map:' mfilename ':mapdispError'], ...
                'Incorrect %s property', upper(propname))
        end

    case 'geoid'
        if ischar(propvalu)
            if exist(propvalu,'file') == 2
                mstruct.geoid = feval(propvalu,'geoid');
            else
                error(['map:' mfilename ':mapdispError'], ...
                    'Incorrect %s property', upper(propname))
            end
        else
            mstruct.geoid = geoidtst(propvalu);
            if mstruct.geoid(1) == 0
                error(['map:' mfilename ':mapdispError'], ...
                    'Positive Geoid radius required')
            end
        end

    case 'mapparallels'
        if ischar(propvalu) || length(propvalu) > 2
            error(['map:' mfilename ':mapdispError'], ...
                'Incorrect %s property', upper(propname))
        elseif mstruct.nparallels == 0
            error(['map:' mfilename ':mapdispError'], ...
                'PARALLELS can not be specified for this projection')
        elseif length(propvalu) > mstruct.nparallels
            error(['map:' mfilename ':mapdispError'], ...
                'Too many PARALLELS for this projection')
        else
            mstruct.mapparallels = propvalu;
        end

    case 'scalefactor'
        if ~isnumeric(propvalu) || length(propvalu) > 1 ||  propvalu == 0
            error(['map:' mfilename ':mapdispError'], ...
                'Incorrect %s property', upper(propname))
        else
            if strcmp(mstruct.mapprojection,'utm') || ...
                    strcmp(mstruct.mapprojection,'ups')
                error(['map:' mfilename ':mapdispError'], ...
                    'SCALEFACTOR cannot be specified for this projection')
            else
                mstruct.scalefactor = propvalu;
            end
        end

    case 'falseeasting'
        if ~isnumeric(propvalu) || length(propvalu) > 1
            error(['map:' mfilename ':mapdispError'], ...
                'Incorrect %s property', upper(propname))
        else
            if strcmp(mstruct.mapprojection,'utm') || ...
                    strcmp(mstruct.mapprojection,'ups')
                error(['map:' mfilename ':mapdispError'], ...
                    'FALSEEASTING cannot be specified for this projection')
            else
                mstruct.falseeasting = propvalu;
            end
        end

    case 'falsenorthing'
        if ~isnumeric(propvalu) || length(propvalu) > 1
            error(['map:' mfilename ':mapdispError'], ...
                'Incorrect %s property', upper(propname))
        else
            if strcmp(mstruct.mapprojection,'utm') || ...
                    strcmp(mstruct.mapprojection,'ups')
                error(['map:' mfilename ':mapdispError'], ...
                    'FALSENORTHING cannot be specified for this projection')
            else
                mstruct.falsenorthing = propvalu;
            end
        end


    %******************
    %  Frame Properties
    %******************

    case 'frame'
        validparm = {'on','off','reset'};
        indx = strmatch(propvalu,validparm);
        if length(indx) == 1
            if indx == 3;   indx = 1;  end  %  Reset becomes on
            mstruct.frame = validparm{indx};
        else
            error(['map:' mfilename ':mapdispError'], ...
                'Incorrect %s property', upper(propname))
        end

    case 'fedgecolor'
        if ischar(propvalu) || ...
                (length(propvalu) == 3 && all(propvalu <= 1 & propvalu >= 0))
            mstruct.fedgecolor = propvalu;
        else
            error(['map:' mfilename ':mapdispError'], ...
                'Incorrect %s property', upper(propname))
        end

    case 'ffacecolor'
        if ischar(propvalu) || ...
                (length(propvalu) == 3 && all(propvalu <= 1 & propvalu >= 0))
            mstruct.ffacecolor = propvalu;
        else
            error(['map:' mfilename ':mapdispError'], ...
                'Incorrect %s property', upper(propname))
        end

    case 'ffill'
        if ~ischar(propvalu)
            mstruct.ffill = max([propvalu,2]);
        else
            error(['map:' mfilename ':mapdispError'], ...
                'Incorrect %s property', upper(propname))
        end

    case 'flinewidth'
        if ~ischar(propvalu)
            mstruct.flinewidth = max([propvalu(:),0]);
        else
            error(['map:' mfilename ':mapdispError'], ...
                'Incorrect %s property', upper(propname))
        end

    %*************************
    %  General Grid Properties
    %*************************

    case 'grid'
        validparm = {'on','off','reset'};
        indx = strmatch(propvalu,validparm);
        if length(indx) == 1
            if indx == 3;   indx = 1;  end  %  Reset becomes on
            mstruct.grid = validparm{indx};
        else
            error(['map:' mfilename ':mapdispError'], ...
                'Incorrect %s property', upper(propname))
        end

    case 'galtitude'
        if ~ischar(propvalu)
            mstruct.galtitude = propvalu(1);
        else
            error(['map:' mfilename ':mapdispError'], ...
                'Incorrect %s property', upper(propname))
        end

    case 'gcolor'
        if ischar(propvalu) || ...
                (length(propvalu) == 3 && all(propvalu <= 1 & propvalu >= 0))
            mstruct.gcolor = propvalu;
        else
            error(['map:' mfilename ':mapdispError'], ...
                'Incorrect %s property', upper(propname))
        end

    case 'glinestyle'
        [lstyle,~,~,err] = colstyle(propvalu);
        error(err) %#ok<ERTAG>
        mstruct.glinestyle = lstyle;
        if isempty(lstyle)
            warning('map:axesm:missingGridLineStyle', ...
                'Missing grid line style.')
        end

    case 'glinewidth'
        if ~ischar(propvalu)
            mstruct.glinewidth = max([propvalu(:),0]);
        else
            error(['map:' mfilename ':mapdispError'], ...
                'Incorrect %s property', upper(propname))
        end


    %**************************
    %  Meridian Grid Properties
    %**************************

    case 'mlineexception'
        if ischar(propvalu)
            error(['map:' mfilename ':mapdispError'], ...
                'Incorrect %s property', upper(propname))
        else
            mstruct.mlineexception = propvalu;
        end

    case 'mlinefill'
        if ~ischar(propvalu)
            mstruct.mlinefill = max([propvalu, 2]);
        else
            error(['map:' mfilename ':mapdispError'], ...
                'Incorrect %s property', upper(propname))
        end

    case 'mlinelimit'
        if ischar(propvalu) || length(propvalu) ~= 2
            error(['map:' mfilename ':mapdispError'], ...
                'Incorrect %s property', upper(propname))
        else
            mstruct.mlinelimit = propvalu;
        end

    case 'mlinelocation'
        if ischar(propvalu)
            error(['map:' mfilename ':mapdispError'], ...
                'Incorrect %s property', upper(propname))
        elseif length(propvalu) == 1
            mstruct.mlinelocation = abs(propvalu);
        else
            mstruct.mlinelocation = propvalu;
        end

    case 'mlinevisible'
        validparm = {'on','off'};
        indx = strmatch(propvalu,validparm);
        if length(indx) == 1
            mstruct.mlinevisible = validparm{indx};
        else
            error(['map:' mfilename ':mapdispError'], ...
                'Incorrect %s property', upper(propname))
        end


    %**************************
    %  Parallel Grid Properties
    %**************************

    case 'plineexception'
        if ischar(propvalu)
            error(['map:' mfilename ':mapdispError'], ...
                'Incorrect %s property', upper(propname))
        else
            mstruct.plineexception = propvalu;
        end

    case 'plinefill'
        if ~ischar(propvalu)
            mstruct.plinefill = max([propvalu, 2]);
        else
            error(['map:' mfilename ':mapdispError'], ...
                'Incorrect %s property', upper(propname))
        end

    case 'plinelimit'
        if ischar(propvalu) || length(propvalu) ~= 2
            error(['map:' mfilename ':mapdispError'], ...
                'Incorrect %s property', upper(propname))
        else
            mstruct.plinelimit = propvalu;
        end

    case 'plinelocation'
        if ischar(propvalu)
            error(['map:' mfilename ':mapdispError'], ...
                'Incorrect %s property', upper(propname))
        elseif length(propvalu) == 1
            mstruct.plinelocation = abs(propvalu);
        else
            mstruct.plinelocation = propvalu;
        end

    case 'plinevisible'
        validparm = {'on','off'};
        indx = strmatch(propvalu,validparm);
        if length(indx) == 1
            mstruct.plinevisible = validparm{indx};
        else
            error(['map:' mfilename ':mapdispError'], ...
                'Incorrect %s property', upper(propname))
        end


    %**************************
    %  General Label Properties
    %**************************

    case 'fontangle'
        validparm = {'normal','italic','oblique'};
        indx = strmatch(propvalu,validparm);
        if length(indx) == 1
            mstruct.fontangle = validparm{indx};
        else
            error(['map:' mfilename ':mapdispError'], ...
                'Incorrect %s property', upper(propname))
        end

    case 'fontcolor'
        if ischar(propvalu) || ...
                (length(propvalu) == 3 && all(propvalu <= 1 & propvalu >= 0))
            mstruct.fontcolor = propvalu;
        else
            error(['map:' mfilename ':mapdispError'], ...
                'Incorrect %s property', upper(propname))
        end

    case 'fontname'
        if ischar(propvalu)
            mstruct.fontname = propvalu;
        else
            error(['map:' mfilename ':mapdispError'], ...
                'Incorrect %s property', upper(propname))
        end

    case 'fontsize'
        if ischar(propvalu) || length(propvalu) ~= 1
            error(['map:' mfilename ':mapdispError'], ...
                'Incorrect %s property', upper(propname))
        else
            mstruct.fontsize = propvalu;
        end

    case 'fontunits'
        validparm = {'points','normalized','inches',...
            'centimeters','pixels'};
        indx = strmatch(propvalu,validparm);
        if length(indx) == 1
            mstruct.fontunits = validparm{indx};
        else
            error(['map:' mfilename ':mapdispError'], ...
                'Incorrect %s property', upper(propname))
        end

    case 'fontweight'
        validparm = {'normal','bold'};
        indx = strmatch(propvalu,validparm);
        if length(indx) == 1
            mstruct.fontweight = validparm{indx};
        else
            error(['map:' mfilename ':mapdispError'], ...
                'Incorrect %s property', upper(propname))
        end

    case 'labelformat'
        validparm = {'compass','signed','none'};
        indx = strmatch(propvalu,validparm);
        if length(indx) == 1
            mstruct.labelformat = validparm{indx};
        else
            error(['map:' mfilename ':mapdispError'], ...
                'Incorrect %s property', upper(propname))
        end

    case 'labelunits'
        if strncmpi(propvalu, 'dms', numel(propvalu))
            mstruct.labelunits = lower(propvalu);
        else
            mstruct.labelunits = checkangleunits(propvalu);
        end

    case 'labelrotation'
        validparm = {'on','off'};
        indx = strmatch(propvalu,validparm);
        if length(indx) == 1
            mstruct.labelrotation = validparm{indx};
        else
            error(['map:' mfilename ':mapdispError'], ...
                'Incorrect %s property', upper(propname))
        end


    %***************************
    %  Meridian Label Properties
    %***************************

    case 'meridianlabel'
        validparm = {'on','off','reset'};
        indx = strmatch(propvalu,validparm);
        if length(indx) == 1
            if indx == 3;   indx = 1;  end  %  Reset becomes on
            mstruct.meridianlabel = validparm{indx};
        else
            error(['map:' mfilename ':mapdispError'], ...
                'Incorrect %s property', upper(propname))
        end

    case 'mlabellocation'
        if ischar(propvalu)
            error(['map:' mfilename ':mapdispError'], ...
                'Incorrect %s property', upper(propname))
        elseif length(propvalu) == 1
            mstruct.mlabellocation = abs(propvalu);
        else
            mstruct.mlabellocation = propvalu;
        end

    case 'mlabelparallel'
        if ischar(propvalu)
            validparm = {'north','south','equator'};
            indx = strmatch(propvalu,validparm);
            if length(indx) == 1
                mstruct.mlabelparallel = validparm{indx};
            else
                error(['map:' mfilename ':mapdispError'], ...
                    'Incorrect %s property', upper(propname))
            end
        elseif length(propvalu) == 1
            mstruct.mlabelparallel = propvalu;
        else
            error(['map:' mfilename ':mapdispError'], ...
                'Incorrect %s property', upper(propname))
        end

    case 'mlabelround'
        if ischar(propvalu) || length(propvalu) ~= 1
            error(['map:' mfilename ':mapdispError'], ...
                'Incorrect %s property', upper(propname))
        else
            mstruct.mlabelround = round(propvalu);
        end


    %***************************
    %  Parallel Label Properties
    %***************************

    case 'parallellabel'
        validparm = {'on','off','reset'};
        indx = strmatch(propvalu,validparm);
        if length(indx) == 1
            if indx == 3;   indx = 1;  end  %  Reset becomes on
            mstruct.parallellabel = validparm{indx};
        else
            error(['map:' mfilename ':mapdispError'], ...
                'Incorrect %s property', upper(propname))
        end

    case 'plabellocation'
        if ischar(propvalu)
            error(['map:' mfilename ':mapdispError'], ...
                'Incorrect %s property', upper(propname))
        elseif length(propvalu) == 1
            mstruct.plabellocation = abs(propvalu);
        else
            mstruct.plabellocation = propvalu;
        end

    case 'plabelmeridian'
        if ischar(propvalu)
            validparm = {'east','west','prime'};
            indx = strmatch(propvalu,validparm);
            if length(indx) == 1
                mstruct.plabelmeridian = validparm{indx};
            else
                error(['map:' mfilename ':mapdispError'], ...
                    'Incorrect %s property', upper(propname))
            end
        elseif length(propvalu) == 1
            mstruct.plabelmeridian = propvalu;
        else
            error(['map:' mfilename ':mapdispError'], ...
                'Incorrect %s property', upper(propname))
        end

    case 'plabelround'
        if ischar(propvalu) || length(propvalu) ~= 1
            error(['map:' mfilename ':mapdispError'], ...
                'Incorrect %s property', upper(propname))
        else
            mstruct.plabelround = round(propvalu);
        end

    otherwise
        error(['map:' mfilename ':mapdispError'], ...
            'Read only property %s can not be modified', upper(propname))
end

%-----------------------------------------------------------------------

function mstruct = setflatlimit(mstruct, flatlimit)

if ischar(flatlimit) || length(flatlimit) > 2
    error(['map:' mfilename ':mapdispError'], ...
        'Incorrect %s property', upper(propname))
elseif strcmp(mstruct.mapprojection,'globe')
    warning('map:axesm:flatlimitGlobe', ...
         ['Ignoring value of %s. %s cannot be set for the ' ...
          '%s projection.'], ...
          'FLatLimit', 'FLatLimit','''globe''')
else
    mstruct.flatlimit = flatlimit;
end

%-----------------------------------------------------------------------

function mstruct = setflonlimit(mstruct, flonlimit)

if ischar(flonlimit) || (length(flonlimit) ~= 2 && ~isempty(flonlimit))
    error(['map:' mfilename ':mapdispError'], ...
        'Incorrect %s property', upper(propname))
elseif strcmp(mstruct.mapprojection,'globe')
    warning('map:axesm:flonlimitGlobe', ...
         ['Ignoring value of %s. %s cannot be set for the ' ...
          '%s projection.'], ...
          'FLonLimit', 'FLonLimit','''globe''')
else
    mstruct.flonlimit = flonlimit;
end

%-----------------------------------------------------------------------

function setgca(mstruct)

%  Set GCA to be map axes.
set(gca, ...
    'NextPlot','Add',...
    'UserData',mstruct,...
    'DataAspectRatio',[1 1 1],...
    'Box','on', ...
    'ButtonDownFcn','uimaptbx')

%  Show the axes background but not the axes labels.
showaxes('off');

%-----------------------------------------------------------------------

function setframegrat(mstruct)

%  Display grid and frame if necessary
if strcmp(mstruct.frame,'on')
    framem('reset');
end

if strcmp(mstruct.grid,'on')
    gridm('reset');
end

if strcmp(mstruct.meridianlabel,'on')
    mlabel('reset');
end

if strcmp(mstruct.parallellabel,'on')
    plabel('reset');
end
