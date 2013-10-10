function [hndl,msg] = handlem(object,axishndl,method)
%HANDLEM Handles of displayed map objects
%
%  HANDLEM or HANDLEM TAGLIST displays a dialog box for selecting objects
%  that have their Tag property set.
%
%  HANDLEM PROMPT displays a dialog box for selecting objects based on the
%  object strings listed below.
%
%  H = HANDLEM(OBJECT) returns the handles of those objects in the current
%  axes specified by the input string, OBJECT. The options for the object
%  string are defined by the following list:
%
%      ALL         All children 
%      CLABEL      Contour labels
%      CONTOUR     hggroups containing contours
%      FILLCONTOUR hggroups containing filled contours
%      FRAME       Map frame
%      GRID        Map grid lines
%      HGGROUP     All hggroup objects
%      HIDDEN      Hidden objects
%      IMAGE       Untagged image objects
%      LIGHT       Untagged light objects
%      LINE        Untagged line objects
%      MAP         All objects on the map, excluding the frame and grid
%      MERIDIAN    Longitude grid lines
%      MLABEL      Longitude labels
%      PARALLEL    Latitude grid lines
%      PLABEL      Latitude labels
%      PATCH       Untagged patch objects
%      SCALERULER  Scaleruler objects
%      SURFACE     Untagged surface objects
%      TEXT        Untagged text objects
%      TISSOT      Tissot indicatrices
%      VISIBLE     Visible objects
%
%  H = HANDLEM(TAGSTR) returns the handles for any objects whose tags match
%  the string TAGSTR.
%
%  H = HANDLEM(OBJECT, AXESH) or HANDLEM(TAGSTR, AXESH) searches within the
%  axes specified by the input handle AXESH.
%
%  H = HANDLEM(..., AXESH, SEARCHMETHOD) controls the method used to match
%  the OBJECT input. If omitted, 'exact' is assumed. Search method
%  'strmatch' searches for matches that start at the beginning of the tag.
%  Search method 'findstr' searches anywhere within the tag for the object
%  string.
%
%  H = HANDLEM(HANDLES) returns those elements of an input vector of
%  handles that are still valid.
%
%  A prefix of 'all' may be applied to strings defining a Handle Graphics
%  object type (text, line, patch, light, surface, or image) to find  all
%  object handles that meet the type criteria (for example, 'allimage').
%  Without the 'all' prefix, only handles with an empty tag are returned.
%
%  See also FINDOBJ.

% Copyright 1996-2010 The MathWorks, Inc.
% $Revision: 1.11.4.9 $  $Date: 2010/06/26 04:57:27 $

% Obsolete syntax
% ---------------
% [h,msg] = HANDLEM(...) returns a string indicating any error encountered.
if nargout > 1
    warnObsoleteMSGSyntax(mfilename)
    msg = '';
end

%  Define recognized names.  Note that this definition process is
%  significantly faster than strvcat, where the padding must be computed
names = [
    'all        '
    'clabel     '
    'contour    '
    'fillcontour'
    'frame      '
    'grid       '
    'hggroup    '
    'hidden     '
    'image      '
    'light      '
    'line       '
    'map        '
    'meridian   '
    'mlabel     '
    'parallel   '
    'patch      '
    'plabel     '
    'scaleruler '
    'surface    '
    'text       '
    'tissot     '
    'visible    '
    ];

%  Initialize input arguments if necessary
if nargin == 0 || isempty(object) || strcmpi(object,'taglist')
    hndl = PromptFromTags;
    return
elseif strcmpi(object,'prompt')
    object = PromptForName(names);
    if isempty(object)
        return
    end
end

if nargin < 3
    method = 'exact';
else
    assert(any(strcmpi(method, {'exact','findstr','strmatch'})), ...
        'map:handlem:mapdispError', ...
        'Search method must be ''exact'',''findstr'', or ''strmatch''.');
end

%  Test if an axis handle has been provided
if nargin < 2
    axishndl = gca;
else
    if ~(length(axishndl) == 1 && ishghandle(axishndl,'axes') )
        error(['map:' mfilename ':mapdispError'], ...
            'Valid axes handle required.')
    end
end

%  Test for valid handles if input is a numeric vector
if ~ischar(object)
    indx = find(ishghandle(object));
    if isempty(indx)
        error(['map:' mfilename ':mapdispError'], ...
            'Object not found on current axes.')
    else
        hndl = object(indx);
    end
    return
end

%  Test for a valid string input
if ~ischar(object) || min(size(object)) ~= 1
    error(['map:' mfilename ':mapdispError'], ...
        'String vector required')
else
    object = object(:)';    
    %  Enforce row string vector
    %  Allow object to retain letter case for
    %  otherwise matches of tags.
    
    %  Test if prefix of all is applied  
    allflag = 0;
    if length(object) >= 3
        if strcmpi(object(1:3),'all') &&  ...
                ~isempty(strmatch(object(4:length(object)),...
                {'image', 'line', 'surface', 'patch', 'text'},'exact'))
            allflag = 1;    object(1:3) = [];
        end
    end
    
    % test for keyword object. Require exact match
    strindx = strmatch(lower(object),names,'exact');      
    keywordmatch = 0; 
    if isempty(strindx)  || ~strcmp(method,'exact')
        name = object;
    elseif length(strindx) == 1 && strcmp(method,'exact')
        name = deblank(names(strindx,:));   %  Set the name string
        keywordmatch = 1;
    elseif length(strindx) > 1
        msg = ['Object not found on current axes:  ',object];
        error(['map:' mfilename ':mapdispError'], msg)
    end
end

%  Get the children of the current axes
children = get(axishndl,'Children');

%  Set the appropriate handle vector. Use an obfuscated method of enforcing
%  new match method to avoid restructuring code.
if keywordmatch   
    switch name
        case 'all'
            hndl = children;
            
        case 'clabel'
            g = findContourHandle(children);
            hndl = findobj(g, 'Type', 'text');
            
        case 'contour'
            hndl = findContourHandle(children);
            
        case 'fillcontour'
            hndl = findFillContourHandle(children);
            
        case 'hggroup'
            hndl = findobj(children,'Type','hggroup');
            
        case 'frame'
            hndl = findobj(children,'Tag','Frame');
            
        case 'grid'
            lathndl  = findobj(children,'Tag','Parallel');
            lonhndl  = findobj(children,'Tag','Meridian');
            hndl = [lathndl; lonhndl];
            if isempty(hndl)
                hndl = findobj(children,'Tag',object);   
            end
            
        case 'hidden'
            hndl = findobj(children,'Visible','off');
            
        case 'image'
            if allflag
                hndl = findobj(children,'Type','image');
            else
                hndl = findobj(children,'Type','image','Tag','');
            end
            
        case 'light'
            if allflag
                hndl = findobj(children,'Type','light');
            else
                hndl = findobj(children,'Type','light','Tag','');
            end
            
        case 'line'
            if allflag
                hndl = findobj(children,'Type','line');
            else
                hndl = findobj(children,'Type','line','Tag','');
            end
            
        case 'map'
            gcm;   %  Enforce a map axes here
            border = findobj(children,'Tag','Frame');
            lathndl  = findobj(children,'Tag','Parallel');
            lonhndl  = findobj(children,'Tag','Meridian');
            hndl = children;
            
            if ~isempty(border);   hndl( hndl == border )  = [];  end
            if ~isempty(lathndl);  hndl( hndl == lathndl ) = [];  end
            if ~isempty(lonhndl);  hndl( hndl == lonhndl ) = [];  end
            
        case 'meridian'
            hndl = findobj(children,'Tag','Meridian');
            
        case 'mlabel'
            hndl = findobj(children,'Tag','MLabel');
            
        case 'patch'
            if allflag
                hndl = findobj(children,'Type','patch');
            else
                hndl = findobj(children,'Type','patch','Tag','');
            end
            
        case 'parallel'
            hndl = findobj(children,'Tag','Parallel');
            
        case 'plabel'
            hndl = findobj(children,'Tag','PLabel');
            
        case 'surface'
            if allflag
                hndl = findobj(children,'Type','surface');
            else
                hndl = findobj(children,'Type','surface','Tag','');
            end
            
        case 'text'
            if allflag
                hndl = findobj(children,'Type','text');
            else
                hndl = findobj(children,'Type','text','Tag','');
            end
            
        case 'tissot'
            hndl = findobj(children,'Tag','Tissot');
            
        case 'visible'
            hndl = findobj(children,'Visible','on');
            
        case 'scaleruler'
            hndl = findall(children, 'Type', 'hggroup', ...
                '-regexp', 'Tag',  'scaleruler*');     
    end
else
    % not keyword match
    namecopy = deblank(name);
    if isempty(namecopy)  
        namecopy = name;  
    end     %  Using tag lists, name may be padded with trailing spaces
    
    switch method        
        case 'exact'
            hndl = findobj(children,'Tag',namecopy) ;
            hggroupHndl = findall(hndl, 'type','hggroup');
            if ~isempty(hggroupHndl)
                tagstr = get(hggroupHndl,'Tag');
                if strmatch('scaleruler',tagstr)
                    hndl = hggroupHndl;
                end
            end
                      
        otherwise            
            % get the associated tags
            tagcellarray = get(children,'tag');
            tagstrmat = char(tagcellarray{:});
            
            % search for matches
            indx = findstrmat(tagstrmat,namecopy,method);
            hndl = children(indx);
    end
    if isempty(hndl)
        msg = ['Unrecognized object name:  ',name];
        error(['map:' mfilename ':mapdispError'], msg)
    end
end

%--------------------------------------------------------------------------

function str = PromptForName(namelist)
%  PROMPTFORNAME will produce a modal dialog box which
%  allows the user to select from the recognized HANDLEM.
%  name list or enter their own object tag name.

str = [];
children = get(gca,'Children');

%  Eliminate objects from list if their handles are hidden.
if length(children) == 1
    if strcmp(get(children,'HandleVisibility'),'callback')
        indx = 1;
    else
        indx = [];
    end
else
    hidden = get(children,'HandleVisibility');
    indx = strmatch('callback',char(hidden));
end
children(indx) = [];   

if isempty(children)
    uiwait(errordlg('No objects on current axes',...
        'Object Specification','modal'));  return
end

objnames = namem(children);   %  Names of objects on current axes

%  Compute the Pixel and Font Scaling Factors so
%  GUI figure windows and fonts look OK across all platforms
PixelFactor = guifactm('pixels');
FontScaling =  guifactm('fonts');

%  Create the dialog box.  Make visible when all objects are drawn
h = dialog('Name','Specify Object',...
    'Units','Points',  'Position',PixelFactor*72*[2 1 3 3.5],...
    'Visible','off');
colordef(h,'white');
figclr = get(h,'Color');

%  Object Name Title and Frame
uicontrol(h,'Style','Frame',...
    'Units','Points',  'Position',PixelFactor*72*[0.03  1.85  2.94  1.45], ...
    'ForegroundColor', 'black','BackgroundColor', figclr);

uicontrol(h,'Style','Text','String','Object', ...
    'Units','Points',  'Position',PixelFactor*72*[0.09  3.26  0.90  0.20], ...
    'FontWeight','bold',  'FontSize',FontScaling*12, ...
    'HorizontalAlignment', 'center', ...
    'ForegroundColor', 'black','BackgroundColor', figclr);

%  Object Name Text and Popup Menu
uicontrol(h,'Style','Text','String', 'Name:', ...
    'Units','Points',  'Position',PixelFactor*72*[0.18  2.90  0.90  0.24], ...
    'FontWeight','bold',  'FontSize',FontScaling*10, ...
    'HorizontalAlignment', 'right', ...
    'ForegroundColor', 'black','BackgroundColor', figclr);

p = uicontrol(h,'Style','Popup','String', namelist,'Value', 1, ...
    'Units','Points',  'Position',PixelFactor*72*[1.20  2.86  1.50  0.32], ...
    'FontWeight','bold',  'FontSize',FontScaling*10, ...
    'HorizontalAlignment', 'center',...
    'ForegroundColor', 'black','BackgroundColor', figclr);

%  Other Tag and Edit Box
uicontrol(h,'Style','Text','String', 'Other Tag:', ...
    'Units','Points',  'Position',PixelFactor*72*[0.18  2.46  0.90  0.24], ...
    'FontWeight','bold',  'FontSize',FontScaling*10, ...
    'HorizontalAlignment', 'right', ...
    'ForegroundColor', 'black','BackgroundColor', figclr);

callbackstr = ['if ~isempty(get(gco,''String''));'                   ,...
    '     set(get(gco,''UserData''),''Enable'',''off'');' ,...
    'else;  set(get(gco,''UserData''),''Enable'',''on'');end'];

e = uicontrol(h,'Style','Edit','String', '', ...
    'Units','Points',  'Position',PixelFactor*72*[1.20  2.42  1.50  0.32], ...
    'FontWeight','bold',  'FontSize',FontScaling*10, ...
    'HorizontalAlignment', 'left', ...
    'ForegroundColor', 'black','BackgroundColor', figclr,...
    'UserData',p,'CallBack',callbackstr);

%  Namelist select button
callbackstr = [
    'get(gco,''UserData'');',...
    'ans.indx = listdlg(''ListString'',cellstr(ans.objnames),''SelectionMode'',''single'',',...
    '''ListSize'',[160 170],''Name'',''Select Object'');',...
    'if ~isempty(ans.indx);  ans.ud = get(gco,''UserData'');',...
    'set(ans.ud.hndl(1),''String'',deblank(ans.objnames(ans.indx,:)));',...
    'set(ans.ud.hndl(2),''Enable'',''off'');end;clear ans'];

userdata.hndl = [e p];   userdata.objnames = objnames;

uicontrol(h,'Style','Push','String', 'Select', ...
    'Units','Points',  'Position',PixelFactor*72*[1.40  2.00  1.10  0.30], ...
    'FontWeight','bold',  'FontSize',FontScaling*10, ...
    'HorizontalAlignment', 'center',...
    'UserData', userdata,...
    'ForegroundColor', 'black', 'BackgroundColor', figclr,...
    'Interruptible','on','CallBack',callbackstr);

%  Match Title and Frame
uicontrol(h,'Style','Frame',...
    'Units','Points',  'Position',PixelFactor*72*[0.03  0.70  2.94  0.90], ...
    'ForegroundColor', 'black','BackgroundColor', figclr);

uicontrol(h,'Style','Text', 'String','Match', ...
    'Units','Points',  'Position',PixelFactor*72*[0.09  1.5  0.90  0.20], ...
    'FontWeight','bold',  'FontSize',FontScaling*12, ...
    'HorizontalAlignment', 'center', ...
    'ForegroundColor', 'black','BackgroundColor', figclr);

%  Match Radio Buttons
r1 = uicontrol(h,'Style','Radio','String','Untagged Objects', 'Value',1,...
    'Units','Points',  'Position',PixelFactor*72*[0.18 1.14 2.00 0.32],...
    'FontWeight','bold',  'FontSize',FontScaling*10, ...
    'HorizontalAlignment', 'left', ...
    'ForegroundColor', 'black','BackgroundColor', figclr);

r2 = uicontrol(h,'Style','Radio', 'String','All Objects', 'Value',0,...
    'Units','Points',  'Position',PixelFactor*72*[0.18 0.85 2.00 0.32],...
    'FontWeight','bold',  'FontSize',FontScaling*10, ...
    'HorizontalAlignment','left', 'Tag','AllObjects', ...
    'ForegroundColor', 'black','BackgroundColor', figclr,...
    'Interruptible','on','Enable','off');

%  Set the user data properties and callbacks for the radio buttons
callbackstr = 'set(gco,''Value'',1);set(get(gco,''UserData''),''Value'',0);';
set(r1,'UserData',r2,'CallBack',callbackstr)
set(r2,'UserData',r1,'CallBack',callbackstr)

%  Buttons to exit the modal dialog
uicontrol(h,'Style','Push','String', 'Apply', ...    %  Apply Button
    'Units','Points',  'Position',PixelFactor*72*[0.30  0.10  1.05  0.40], ...
    'FontWeight','bold',  'FontSize',FontScaling*12, ...
    'HorizontalAlignment', 'center', ...
    'ForegroundColor', 'black', 'BackgroundColor', figclr,...
    'CallBack','uiresume');

uicontrol(h,'Style','Push','String', 'Cancel', ...    %  Cancel Button
    'Units','Points',  'Position',PixelFactor*72*[1.65  0.10  1.05  0.40], ...
    'FontWeight','bold',  'FontSize',FontScaling*12, ...
    'HorizontalAlignment', 'center', ...
    'ForegroundColor', 'black','BackgroundColor', figclr,...
    'CallBack','uiresume');

%  Set the callback for the popup menu.  Disable all match option
%  if not a handle graphics child of an axes.
callbackstr = [
    'get(gco,''String'');deblank(ans(get(gco,''Value''),:));',...
    'if strcmp(ans,''text'') | strcmp(ans,''patch'') | '     ,...
    'strcmp(ans,''surface'') | strcmp(ans,''image'') | '     ,...
    'strcmp(ans,''light'') | strcmp(ans,''line'');'          ,...
    'get(gco,''UserData'');set(ans(2),''Enable'',''on'');'   ,...
    'else;get(gco,''UserData'');set(ans(1),''Value'',1);'    ,...
    'set(ans(2),''Enable'',''off'',''Value'',0);end;clear ans'];

set(p,'UserData',[r1 r2],'CallBack',callbackstr)

%  Turn dialog box on.  Then wait unit a button is pushed
set(h,'Visible','on');     uiwait(h)

if ~ishghandle(h)
    return;
end

%  If the accept button has been pushed, then
%  first determine if the object edit box has a string in
%  it.  If it does not, then get the name from the
%  popup menu with the name list.  Finally, check
%  to see if the all match option is selected.  If so,
%  append 'all' to the string.
if strcmp(get(get(h,'CurrentObject'),'String'),'Apply')
    str = get(e,'String');
    if isempty(str)
        str = deblank(namelist(get(p,'Value'),:));   
    end
    if get(r2,'Value')
        str = ['All', str];   
    end
end

%  Close the dialog box
delete(h)

%--------------------------------------------------------------------------

function hndl = PromptFromTags
%  PROMPTFROMTAGS will produce a modal dialog box which
%  allows the user to select from the recognized object TAGS on
%  the current axes.

hndl = [];
children = get(gca,'Children');
if isempty(children);
    uiwait(errordlg('No objects on current axes',...
        'Object Specification','modal'));  return;
end

%  Eliminate objects from list if their handles are hidden.
if length(children) == 1
    if strcmp(get(children,'HandleVisibility'),'callback')
        indx = 1;
    else
        indx = [];
    end
else
    hidden = get(children,'HandleVisibility');
    indx = strmatch('callback',char(hidden));
end
if ~isempty(indx);   children(indx) = [];   end

%  Display the list dialog if children remain without hidden handles
if ~isempty(children)
    objnames = namem(children);
    indx = listdlg('ListString',cellstr(objnames),...
        'SelectionMode','multiple',...
        'ListSize',[160 170],...
        'Name','Select Object');
    if ~isempty(indx)
        for i = 1:length(indx)
            hndl0 = handlem(deblank(objnames(indx(i),:)));
            hndl = [hndl;hndl0]; %#ok<AGROW>
        end
    end
    
else
    uiwait(errordlg('No objects on current axes',...
        'Object Specification','modal'));  return;
end

%--------------------------------------------------------------------------

function indx = findstrmat(strmat,searchstr,method)

strmat(:,end+1) = 13; % add a lineending character to prevent matches across rows

% find matches in vector
switch method
    case 'findstr'
        % make string matrix a vector
        sz = size(strmat);
        strmat = strmat';
        strvec = strmat(:)';
        vecindx = findstr(searchstr,strvec);
        % vector indices to row indices
        indx = unique(ceil(vecindx/sz(2)));
    case 'strmatch'
        indx = strmatch(searchstr,strmat);
    case 'exact'
        % added a lineending character above to prevent matches across rows
        searchstr(end+1) = 13; 
        indx = strmatch(searchstr,strmat,'exact'); %#ok<*REMFF1>
end

%--------------------------------------------------------------------------

function hndl = findContourHandle(children)
% Find contour hggroup handles.

hndl = findobj(children, 'Type', 'hggroup');
if ~isempty(hndl)
    tf = false(size(hndl));
    for k=1:numel(tf)
        h = hndl(k);
        if isappdata(h,'mapgraph')
            obj = getappdata(h, 'mapgraph');
            if ~isempty(findprop(obj, 'Fill'))
                tf(k) = true;
            end
        end
    end
    hndl = hndl(tf);
end

%--------------------------------------------------------------------------

function hndl = findFillContourHandle(children)
% Find filled contour hggroup handles.

hndl = findContourHandle(children);
tf = false(size(hndl));
for k=1:numel(tf)
    h = hndl(k);
    obj = getappdata(h, 'mapgraph');
    if isequal(obj.Fill, 'on')
        tf(k) = true;
    end
end
hndl = hndl(tf);
