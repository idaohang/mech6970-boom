function vmap0ui(devicename)
%VMAP0UI UI for selecting data from Vector Map Level 0
%
%   VMAP0UI(DIRNAME) launches a graphical user interface for
%   interactively selecting and importing data from a Vector Map Level 0
%   (VMAP0) data base.  Use the string DIRNAME to specify the directory
%   containing the data base.  For more on using VMAP0UI, click the HELP
%   button after the interface appears.
%
%   VMAP0UI(DEVICENAME) or VMAP0UI DEVICENAME uses the logical device
%   (volume) name specified in string DEVICENAME to locate CD-ROM drive
%   containing the VMAP0 CD-ROM.  Under the Windows operating system it 
%   could be 'F:' or 'G:' or some other letter.  Under Macintosh OS X it
%   should be '/Volumes/VMAP'.  Under other UNIX systems it could be 
%   '/cdrom/'.
%
%   VMAP0UI can be used on Windows without any arguments.  In this case it
%   attempts to automatically detect a drive containing a VMAP0 CD-ROM. If
%   VMAP0UI fails to locate the CD-ROM device, then specify it explicitly.
%
%   Vector Map Level 0, created in the 1990s, is still probably the most
%   detailed global database of vector map data available to the public.
%   VMAP0 CD-ROMs are available from through the U.S. Geological Survey
%   (USGS):
% 
%       USGS Information Services (Map and Book Sales)
%       Box 25286
%       Denver Federal Center
%       Denver, CO 80225
%       Telephone: (303) 202-4700
%       Fax: (303) 202-4693
%
%   Examples
%   --------
%   % Launch VMAP0UI and automatically detect a CD-ROM on Windows
%   vmap0ui
%
%   % Launch VMAP0UI on Macintosh OS X (need to specify volume name)
%   vmap0ui('Volumes/VMAP')
%
%   See also DISPLAYM, EXTRACTM, MLAYERS, VMAP0DATA.

% Copyright 1996-2010 The MathWorks, Inc.
% $Revision: 1.1.6.8 $  $Date: 2010/08/07 07:28:21 $
% Written by:  W. Stumpf

% Try to autodetect CD if no device name provided.
if nargin == 0
    devicename = getdevicename;
    
    if isempty(devicename)
        error('map:vmap0ui:deviceNotFound', ...
            'Couldn''t detect the VMAP0 CD. Please provide a device name.')
    end
    
    if numel(devicename) > 1
        warning('map:vmap0ui:multipleDataCDsFound',...
            'More than one VMAP0 CD found. Using %s', devicename{1})
    end
    
    devicename = devicename{1};
end
vmap0uistart(devicename)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function foundpth = getdevicename
%GETDEVICENAME searches for VMAPLVL0 directory
%
% foundpath = GETDEVICENAME PC device names for the vmaplvl0 directory. 
% Now try all letter drives on PCs (A-Z)

s.name = 'vmaplv0';
foundpth = {};
if ispc
   
   for j=1:26
      
      pth = [char('A'+j-1) ':\'];
      try 
         
         d = dir(pth);
         d(~[d.isdir]) = [];
         for i=1:length(s)
            if ~isempty(d) && ismember(s(i).name,lower({d.name}))
               foundpth{end+1} = pth;
            end
         end
         
      catch
         %oops, that drive didn't exist
      end
      
   end
   
elseif strcmp(computer,'MAC2')
   pth = 'vmap:';
   try 
      
      d = dir(pth);
      d(~[d.isdir]) = [];
      for i=1:length(s)
         if ~isempty(d) && ismember(s(i).name,lower({d.name}))
            foundpth{end+1} = pth;
         end
      end
      
   catch
      %oops, that drive didn't exist
   end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function vmap0uistart(devicename)

% check for valid inputs


% Check that the top of the database file hierarchy is visible,
% and note the case of the directory and filenames.

filepath = fullfile(devicename,filesep);
dirstruc = dir(filepath);

if isempty(strmatch('VMAPLV0',upper(strvcat(dirstruc.name)),'exact'))
   error('map:vmap0ui:deviceNotMounted1', ...
       'VMAP Level 0 disk not mounted or incorrect devicename.'); 
end

if ~isempty(strmatch('VMAPLV0',{dirstruc.name},'exact'))
   filesystemcase = 'upper';
elseif ~isempty(strmatch('vmaplv0',{dirstruc.name},'exact'))
   filesystemcase = 'lower';
else
   error('map:vmap0ui:mixedCaseFilenames', ...
       'Unexpected mixed case filenames.')
end

% Read the library attribute table to get the library associated with the device

switch filesystemcase
case 'upper'
   filepath = fullfile(devicename,'VMAPLV0',filesep);
case 'lower'
   filepath = fullfile(devicename,'vmaplv0',filesep);
end
LAT = vmap0read(filepath,'LAT');

library = LAT(end).library_name;

% build the pathname so that [pathname filename] is the full filename

switch filesystemcase
case 'upper'
   filepath = fullfile(devicename,'VMAPLV0',upper(library),filesep);
case 'lower'
   filepath = fullfile(devicename,'vmaplv0',lower(library),filesep);
end

dirstruc = dir(filepath);
if isempty(dirstruc)
   error('map:vmap0ui:deviceNotMounted2', ...
       'VMAP Level 0 disk %s not mounted or incorrect devicename.', ...
       upper(library)); 
end

% Read the list of themes on the CD

CAT = vmap0read(filepath,'CAT');

% Remove the metadata themes from the Coverage Attributes Table

indx = find(ismember({CAT.coverage_name},{'libref','tileref','dq'}));
CAT(indx) = [];

% Get the essential libref information (tile name/number, bounding boxes)

switch filesystemcase
case 'upper'
   filepath = fullfile(devicename,'VMAPLV0',upper(library),'TILEREF',filesep);
case 'lower'
   filepath = fullfile(devicename,'vmaplv0',lower(library),'tileref',filesep);
end

tFT = vmap0read(filepath,'TILEREF.AFT');
FBR = vmap0read(filepath,'FBR');


%get information on type and attributes of data within each theme

for i=1:length(CAT)
   
   switch filesystemcase
   case 'upper'
      filepath = fullfile(devicename,'VMAPLV0',upper(library),upper(CAT(i).coverage_name),filesep);
      ifilename = 'INT.VDT';
      cfilename = 'CHAR.VDT';
   case 'lower'
      filepath = fullfile(devicename,'vmaplv0',lower(library),lower(CAT(i).coverage_name),filesep);
      ifilename = 'int.vdt';
      cfilename = 'char.vdt';
   end
   
   if sign(exist([filepath ifilename],'file'))
      [ivdt{i},ivdtheader{i}] = vmap0read(filepath,ifilename);
   end
   if sign(exist([filepath cfilename],'file'))
      [cvdt{i},cvdtheader{i}] = vmap0read(filepath,cfilename);
   end
   
   wildcard = '*FT';
   if strcmp(filesystemcase,'lower'); wildcard = lower(wildcard); end
   FTfilenames{i} = dir([filepath wildcard]);
   
   for j=1:length(FTfilenames{i})
      headerstr = vmap0rhead(filepath,FTfilenames{i}(j).name);
      [field,varlen,description,narrativefile] =  vmap0phead(headerstr);
      FTfilenames{i}(j).description = strrep(description,' Feature Table','');
      FTfilenames{i}(j).fields = field;
   end
   
end

% construct panel and exit

vmap0uiPanel(CAT,FBR,ivdt,ivdtheader,cvdt,cvdtheader,FTfilenames,devicename,library)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function vmap0uiPanel(CAT,FBR,ivdt,ivdtheader,cvdt,cvdtheader,FTfilenames,devicename,library)
%VMAP0UIPANEL constructs the VMAP0UI panel

h.devicename = devicename;
h.library = library;
h.FBR = FBR;
h.CAT = CAT;
h.ivdt = ivdt;
h.cvdt = cvdt;
h.ivdtheader = ivdtheader;
h.cvdtheader = cvdtheader;
h.FTfilenames = FTfilenames;

for i=1:length(CAT)
   h.VMAPnames{i}=h.CAT(i);
end


% Build the GUI panel

%  Create the dialog window

h.fig = figure('Color',[0.8 0.8 0.8], ...
    'units','character',...
	'Position', [5.2000    2.7692  112.2000   33.1538],...
   'Tag','VMAP0UI','Visible','off',...
   'Menubar','none','Name','VMAP0UI','NumberTitle','off');

colordef(h.fig,'white');
figclr = get(h.fig,'Color');
frameclr = brighten(figclr,0.5);

h.list = uicontrol('Style','list','units','Character',...
   'Position',[2.2440    4.9731   33.6600   25.5285],'Max',1E6,...
   'Interruptible','off','BackgroundColor','w');
h.selected = 1;

h.listlabel = uicontrol('Style','text','Units','Character',...
   'Position',[2.2440   30.8331   33.6600    0.9946],'String','Available Data  ', ...
   'FontWeight','bold','backgroundColor',get(h.fig,'Color'));

strs = str2mat(CAT.description);
h.objparent = zeros(1,length(CAT));
h.objhndl = 1:length(CAT);

suffixstr = ' +   ';
strs = [repmat(suffixstr,length(CAT),1) strs];

h.expandable = ones(1,length(CAT));

set(h.list,'String',strs)
set(h.list,'callback','')
set(h.list,'callback', @vmap0uilist)


% buttons

	h.helpbtn = uicontrol('Parent',h.fig, ...
	'Units','character', ...
	'Position',[4.7206    0.8369   13.2178    2.7898], ...
	'Style','push', ...
	'String','Help', ...
   'Tag','helpbtn',...
   'Callback', @vmap0uihelp,...
	'BackgroundColor',frameclr, 'ForegroundColor','black');

	h.clearbtn = uicontrol('Parent',h.fig, ...
	'Units','character', ...
	'Position', [47.2063    0.8369   13.2178    2.7898], ...
	'Style','push', ...
	'String','Clear', ...
	'Tag','clearbtn',...
    'Callback', @vmap0uiclear,...
	'BackgroundColor',frameclr, 'ForegroundColor','black');

	h.getbtn = uicontrol('Parent',h.fig, ...
	'Units','character', ...
	'Position',[18.8825    0.8369   13.2178    2.7898], ...
	'Style','push', ...
	'String','Get', ...
	'Tag','getbtn',...
   'Callback', @vmap0uiget,...
   'Interruptible','on',...
	'BackgroundColor',frameclr, 'ForegroundColor','black');

	h.savebtn = uicontrol('Parent',h.fig, ...
	'Units','character', ...
	'Position',[75.5300    0.8369   13.2178    2.7898], ...
	'Style','push', ...
	'String','Save', ...
	'Tag','getbtn',...
     'Callback', @vmap0uisave,...
	'BackgroundColor',frameclr, 'ForegroundColor','black');

	h.closebtn = uicontrol('Parent',h.fig, ...
	'Units','character', ...
	'Position',[89.6919    0.8369   13.2178    2.7898], ...
	'Style','push', ...
	'String','Close', ...
    'Callback', @vmap0uiclose,...
	'Tag','getbtn',...
	'BackgroundColor',frameclr, 'ForegroundColor','black');


% Overview map


lonlim = [min([FBR.xmin]) max([FBR.xmax])];
latlim = [min([FBR.ymin]) max([FBR.ymax])];

h.axes = axes('units','Character',...
	'Position',[39.2700    4.9731   67.3200   26.5231],...
	'Tag','samplemapaxes');
    
    

% base map data
axesm('pcarree','MapLatlimit',latlim,'MapLonLimit',lonlim)
tightmap
h.lim = [xlim ylim];

framem;
set(handlem('alltext'),'clipping','on')
load('coast');
hold on
geoshow(lat, long, 'Color',.65*[1 1 1]);
clear lat long

zdatam('allline',-.75)

h.tiles = patchesm([FBR.ymin; FBR.ymax; FBR.ymax; FBR.ymin; FBR.ymin],...
				[FBR.xmin; FBR.xmin; FBR.xmax; FBR.xmax; FBR.xmin],'y',...
            'Tag','tiles','edgecolor',.85*[ 1 1 1]);
zdatam(h.tiles,-1);restack(h.tiles,'bottom')
set(h.tiles,'facecolor','y')

h.baseobjects = get(h.axes,'Children');

zoom

set(gcf,'WindowButtonMotionFcn', @windowButtonMotionCallback)


% display count of tiles covered by current view

do = vmap0do(FBR,latlim,lonlim);
set(get(h.axes,'Title'),'String',[ num2str(length(do)) ' tiles'])

% make figure resizeable

set([...
    h.helpbtn,...
    h.clearbtn,...
    h.getbtn,...
    h.savebtn,...
    h.closebtn,...
    h.list,...
    h.listlabel,...
    h.axes],'units','normalized')
    
set([...
    h.helpbtn,...
    h.clearbtn,...
    h.getbtn,...
    h.savebtn,...
    h.closebtn,...
    h.listlabel,...
    ],'fontunits','normalized')
    

% initialize saved VMAP data
h.extracteddata = [];

set(h.fig,'userdata',h,'Visible','on')
set(handlem('tiles'),'facecolor',[0.99          1.00          0.85])

set(h.fig,'HandleVis','callback')

%------------------------------------------------------------------------

function windowButtonMotionCallback(fig, ~)

if overaxes(fig)
    % Change the cursor to a magnifying glass over the map axes
    set(gcf,'pointer','custom', ...
        'PointerShapeCdata',magcursor(),'PointerShapeHotSpot',[7 7])
else
    set(gcf,'pointer','arrow')
end
    
vmap0uitilecount(fig)


function C = magcursor()

C = [
    0 0 0 0 1 1 1 1 0 0 0 0 0 0 0 0
    0 0 1 1 1 0 0 1 1 1 0 0 0 0 0 0
    0 1 1 0 0 0 0 0 0 1 1 0 0 0 0 0
    0 1 0 0 0 0 0 0 0 0 1 0 0 0 0 0
    1 1 0 0 0 0 0 0 0 0 1 1 0 0 0 0
    1 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0
    1 1 0 0 0 0 0 0 0 0 1 1 0 0 0 0
    0 1 0 0 0 0 0 0 0 1 1 0 0 0 0 0
    0 1 1 0 0 0 0 0 0 1 0 0 0 0 0 0
    0 0 1 1 1 0 0 1 1 1 1 0 0 0 0 0
    0 0 0 0 1 1 1 1 0 1 1 1 0 0 0 0
    0 0 0 0 0 1 0 0 0 0 1 1 1 0 0 0
    0 0 0 0 0 0 0 0 0 0 0 1 1 1 0 0
    0 0 0 0 0 0 0 0 0 0 0 0 1 1 1 0
    0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 1
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1
    ];

C(C == 0)= NaN;


function tf = overaxes(fig)
% If the figure fig contains a single axes (which we expect to be the case),
% return true if the cursor is over that axes.

tf = false;
ax = findobj(fig,'Type','axes');
if isscalar(ax)
    % We expect ax to scalar, but check just to be safe, and do nothing
    % if it does not have exactly one element.
    
    % Get cursor location and figure position in the same units as root.
    p = get(0,'PointerLocation');
    units = get(fig,'Units');
    set(fig,'Units',get(0,'Units'))
    figPos = get(fig,'Position');
    set(fig,'Units',units)
    
    % Normalize cursor location relative to figure.
    x = (p(1) - figPos(1)) / figPos(3);
    y = (p(2) - figPos(2)) / figPos(4);
    
    % Get axes position in normalized units.
    units = get(ax,'Units');
    set(ax,'Units','norm')
    r = get(ax,'Position');
    set(ax,'Units',units)
    
    % Compare normalized positions.
    tf = r(1) <= x && x <= r(1) + r(3) ...
      && r(2) <= y && y <= r(2) + r(4);    
end


function vmap0uitilecount(fig)
% Updates the count of visible tiles

h = get(fig,'UserData');
if ~(all(h.lim == [xlim ylim]))
    [latlim,lonlim] = minvtran(gcm(h.axes),xlim,ylim);
    if lonlim(1) > lonlim(2);
        lonlim(1) = lonlim(1)-360;
    end
    
    do = vmap0do(h.FBR,latlim,lonlim);

    if length(do) == 1
        set(get(h.axes,'Title'),'String',[ num2str(length(do)) ' tile'])

    else
        set(get(h.axes,'Title'),'String',[ num2str(length(do)) ' tiles'])
    end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function vmap0uiclear(~,~)
% Removes extracted data from the plot and storage

%check if a get is in progress
if ~isempty( findall(0,'type','figure','tag','TMWWaitbar','UserData',get(gcbo,'Parent')) );
   return
end

h = get(get(gcbo,'Parent'),'userdata');
hobj = get(h.axes,'Children');
hnew = setdiff(hobj,h.baseobjects);
delete(hnew)
h.extracteddata = [];

zdatam(h.tiles,-1)
set(h.tiles,'Facecolor',[0.99          1.00          0.85],'Edgecolor',.85*[ 1 1 1])

set(h.fig,'userdata',h)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function vmap0uiclose(~,~)
% Closes the VMAP0UI panel

%check if a get is in progress
if ~isempty( findall(0,'type','figure','tag','TMWWaitbar','UserData',get(gcbo,'Parent')) );
   return
end

h = get(get(gcbo,'Parent'),'userdata');
close(gcbf)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function vmap0uisave(~,~)
% Saves VMAP0UI data to a MAT-File

%check if a get is in progress
if ~isempty( findall(0,'type','figure','tag','TMWWaitbar','UserData',get(gcbo,'Parent')) );
   return
end

h = get(get(gcbo,'Parent'),'userdata');

if isempty(h.extracteddata)
   warndlg('No data to save','VMAP0UI warning')
   return
end

answer = questdlg('Save data to a Mat-file or the base workspace?', ...
   'VMAP0UI Save',...
   'Mat-File','Workspace','Cancel','Mat-file');

switch answer
case 'Mat-File'
   vmap0filesave(h.extracteddata)
case 'Workspace'
   vmap0workspacesave(h.extracteddata)
case 'Cancel' 
   return
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function vmap0workspacesave(extracteddata)

answer = questdlg(...
    		strvcat({' ',...
      			'The following variables will be ',...
      			'created in the base workspace:',' ',...
      			extracteddata{:,2},' '}),...
   			'VMAP0UI SAVE','OK','Cancel','OK');

switch answer
case 'OK'
   for i=1:size(extracteddata,1)
      assignin('base',extracteddata{i,2},extracteddata{i,1})
   end
case 'Cancel'
   return
end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function vmap0filesave(VMAPdata)

curpwd = pwd;
[filename,pathname] = uiputfile('*.mat', 'Save the VMAP data in MAT-file:');

if filename == 0
   return
end


eval([ VMAPdata{1,2} ' = VMAPdata{1,1};' ]);
save([pathname filename],VMAPdata{1,2})
clear(VMAPdata{1,2});

for i=2:size(VMAPdata,1)
   
   eval([ VMAPdata{i,2} ' = VMAPdata{i,1};' ]);
   save([pathname filename],VMAPdata{i,2},'-APPEND')
   clear(VMAPdata{i,2});
   
end

cd(curpwd)
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function vmap0uiget(~,~)
% Extracts data from the VMAP0 CD

%check if another get is in progress
if ~isempty( findall(0,'type','figure','tag','TMWWaitbar','UserData',get(gcbo,'Parent')) );
   return
end

% Retrieve handles and such from the GUI user data slot

h = get(get(gcbo,'Parent'),'userdata');

% Determine latitude and longitude limits from current axis limits

x = get(h.axes,'XLim');
y = get(h.axes,'Ylim');
[latlim,lonlim] = minvtran([min(x) max(x)],[min(y) max(y)]);
if lonlim(2) < lonlim(1); lonlim(1) = lonlim(1)-360; end

% Determine which tiles are visible

do = vmap0do(h.FBR,latlim,lonlim);

% Check if the user really wants to extract ALL the tiles

if length(do) == length(h.FBR)-1
   answer = questdlg(['Are you sure you want the entire area covered by the CD?' ...
         ' You can use the mouse to zoom into a smaller region'],...
      ' VMAP0UI Warning','Cancel','Continue','Cancel');
   
   if strcmp(answer,'Cancel')
      return
   end
   
   elseif length(do) > 20
   answer = questdlg(['Are you sure you want to extract ' num2str(length(do)) ' tiles?' ...
         ' You can use the mouse to zoom into a smaller region'],...
      ' VMAP0UI Warning','Cancel','Continue','Cancel');
   
   if strcmp(answer,'Cancel')
      return
   end
   
elseif isempty(do)
   
   warndlg(['You chose an area outside this CD-ROM''s coverage. '...
         'Zoom to another area or close VMAP0UI and try another CD-ROM'],...
         'VMAP0UI Warning');
   
   return
end

% Get pseudohandles for the currently selected lines of the listbox

indx=get(h.list,'value');    									% rows of listbox
objhndl = h.objhndl;												% pseudohandles of each line in the listbox
objparent = h.objparent;										% pseudohandles of the parents of each line in the listbox
hndls = objhndl(indx);											% pseudohandles of selected objects in list
parents = objparent(indx);											% pseudohandles of parents of selected objects in list
VMAPnames = h.VMAPnames(indx);										% contains structures with VMAP file, property and value names

n = gennum(h.objhndl(indx),h.objhndl,h.objparent); 	% generation numbers of the selected rows

% Remove items that are children of selected objects. If the parent is selected, 
% all of the children will also be extracted.

chil = allekinder(hndls,objhndl,objparent);
indx2remove=find(ismember(hndls,chil));
hndls(indx2remove)=[];
parents(indx2remove)=[];
indx(indx2remove)=[];
n(indx2remove)=[];
VMAPnames(indx2remove)=[];

% Assemble a list of themes, feature, properties and values in the
% combinations that would be used with VMAP0DATA.

[themeReq,featureReq,propvalReq,nsteps]=assembleRequests(indx,n,VMAPnames,objhndl,objparent,h);

% get requested themes
hwaitbar = waitbar(0,'Reading data from the VMAP0 CD-ROM...', ...
    'CreateCancelBtn', 'setappdata(gcbf,''canceled'',1)');
setappdata(hwaitbar,'canceled',0)
set(hwaitbar,'handlevisibility','off')
drawnow

nsteps = length(do)*nsteps;
nstepsdone = 0;


try
   for j = 1:length(do)
      
      % obscure the base map data for the current tile
      zdatam(h.tiles(do(j)-1),-.5) 
      set(h.tiles(do(j)-1),'Facecolor','w','EdgeColor','none')
      
      for i=1:length(themeReq)
         
         
         thislatlim = [h.FBR(do(j)).ymin+epsm h.FBR(do(j)).ymax-epsm];
         thislonlim = [h.FBR(do(j)).xmin+epsm h.FBR(do(j)).xmax-epsm];
         
         if isempty(featureReq{i})
            
            % get all features in a requested theme
            [h, nstepsdone] = gettheme(h,hwaitbar,thislatlim,thislonlim,themeReq{i},nstepsdone,nsteps);
            
         elseif isempty(propvalReq{i})
            
            % get just feature table data
            [h, nstepsdone] = getfeature(h,hwaitbar,thislatlim,thislonlim,themeReq{i},featureReq{i},nstepsdone,nsteps);
            
         else
            
            % selected property value pairs of data
            [h, nstepsdone] = getpropval(h,hwaitbar,thislatlim,thislonlim,themeReq{i},featureReq{i},propvalReq{i},nstepsdone,nsteps);
            
         end
      end
   end
catch e
   errordlg(e.message,'Error reading data')
end

if ishghandle(hwaitbar,'figure')
   delete(hwaitbar)
end

% save extracted data back into the figure user data slot

set(h.fig,'userdata',h)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [h,nstepsdone] = getpropval(h,hwaitbar,thislatlim,thislonlim,theme,feature,propval,nstepsdone,nsteps)

feature = lower(feature);


switch feature{:}(end-2:end)
case 'aft'
   topolevel = 'patch';
case 'lft'
   topolevel = 'line';
case 'pft'
   topolevel = 'point';
case 'tft'
   topolevel = 'text';
otherwise
   nstepsdone = nstepsdone+1;
   if ~ishghandle(hwaitbar) || ~ishghandle(gcbo) || getappdata(hwaitbar,'canceled')
      if ishghandle(hwaitbar,'figure'); close(hwaitbar); end
      return;
   else
      waitbar(nstepsdone/nsteps)
      return
   end
end

feature{:}(end-3:end) = [];

% get data:

% Patch

switch topolevel
case 'patch'
   
   VMAPpatch = ...
      vmap0data(h.devicename,h.library,thislatlim,thislonlim,theme,...
      'patch',feature,propval);
   
   if ~ishghandle(hwaitbar) || ~ishghandle(gcbo) || getappdata(hwaitbar,'canceled')
      if ishghandle(hwaitbar,'figure'); close(hwaitbar); end
      return;
   else
      waitbar(nstepsdone/nsteps)
   end
   
   % Set styles of objects to generally coincide with printed OPNAV charts
   [VMAPpatch,newline,newpoint,newtext] = vmap0styles(VMAPpatch);

   h.extracteddata = catmlayerscell(h.extracteddata,VMAPpatch,[feature{:} '_subset']);
   h.extracteddata = catmlayerscell(h.extracteddata,newline,[feature{:} '_subset_line']);
   h.extracteddata = catmlayerscell(h.extracteddata,newpoint,[feature{:} '_subset_patch']);
   h.extracteddata = catmlayerscell(h.extracteddata,newtext,[feature{:} '_subset_text']);
   
   nstepsdone = nstepsdone+1;
   
   if ~ishghandle(hwaitbar) || ~ishghandle(gcbo) || getappdata(hwaitbar,'canceled')
      if ishghandle(hwaitbar,'figure'); close(hwaitbar); end
      return;
   else
      waitbar(nstepsdone/nsteps)
   end
      
   if ~isempty(VMAPpatch); hpatch = displaym(VMAPpatch);zdatam(hpatch,-.25); end
   displaym(newline);
   displaym(newpoint);
   displaym(newtext);

case 'line'
   
   VMAPline = ...
      vmap0data(h.devicename,h.library,thislatlim,thislonlim,theme,...
      'line',feature,propval);
   
   if ~ishghandle(hwaitbar) || ~ishghandle(gcbo) || getappdata(hwaitbar,'canceled')
      if ishghandle(hwaitbar,'figure'); close(hwaitbar); end
      return;
   else
      waitbar(nstepsdone/nsteps)
   end
   
   [VMAPline,newline,newpoint,newtext]= vmap0styles(VMAPline);

   h.extracteddata = catmlayerscell(h.extracteddata,VMAPline,[feature{:} '_subset']);
   h.extracteddata = catmlayerscell(h.extracteddata,newline,[feature{:} '_subset_line']);
   h.extracteddata = catmlayerscell(h.extracteddata,newpoint,[feature{:} '_subset_patch']);
   h.extracteddata = catmlayerscell(h.extracteddata,newtext,[feature{:} '_subset_text']);
   
   nstepsdone = nstepsdone+1;
   
   if ~ishghandle(hwaitbar) || ~ishghandle(gcbo) || getappdata(hwaitbar,'canceled')
      if ishghandle(hwaitbar,'figure'); close(hwaitbar); end
      return;
   else
      waitbar(nstepsdone/nsteps)
   end
   
   if ~isempty(VMAPline); hline = displaym(VMAPline); end
   displaym(newline);
   displaym(newpoint);
   displaym(newtext);

case 'point'
   
   VMAPpoint = ...
      vmap0data(h.devicename,h.library,thislatlim,thislonlim,...
      theme,'point',feature,propval);
   
   if ~ishghandle(hwaitbar) || ~ishghandle(gcbo) || getappdata(hwaitbar,'canceled')
      if ishghandle(hwaitbar,'figure'); close(hwaitbar); end
      return;
   else
      waitbar(nstepsdone/nsteps)
   end
   
   [VMAPpoint,newline,newpoint,newtext]= vmap0styles(VMAPpoint);

   h.extracteddata = catmlayerscell(h.extracteddata,VMAPpoint,[feature{:} '_subset']);
   h.extracteddata = catmlayerscell(h.extracteddata,newline,[feature{:} '_subset_line']);
   h.extracteddata = catmlayerscell(h.extracteddata,newpoint,[feature{:} '_subset_patch']);
   h.extracteddata = catmlayerscell(h.extracteddata,newtext,[feature{:} '_subset_text']);
   
   nstepsdone = nstepsdone+1;
   
   if ~ishghandle(hwaitbar) || ~ishghandle(gcbo) || getappdata(hwaitbar,'canceled')
      if ishghandle(hwaitbar,'figure'); close(hwaitbar); end
      return;
   else
      waitbar(nstepsdone/nsteps)
   end
      
   if ~isempty(VMAPpoint); displaym(VMAPpoint);drawnow; end
   displaym(newline);
   displaym(newpoint);
   displaym(newtext);
   
case 'text'
   
   VMAPtext = ...
      vmap0data(h.devicename,h.library,thislatlim,thislonlim,...
      theme,'text',feature,propval);
   
   if ~ishghandle(hwaitbar) || ~ishghandle(gcbo) || getappdata(hwaitbar,'canceled')
      if ishghandle(hwaitbar,'figure'); close(hwaitbar); end
      return;
   else
      waitbar(nstepsdone/nsteps)
   end
   
   h.extracteddata = catmlayerscell(h.extracteddata,VMAPtext,[feature{:} '_subset']);
   nstepsdone = nstepsdone+1;
   
   if ~ishghandle(hwaitbar) || ~ishghandle(gcbo) || getappdata(hwaitbar,'canceled')
      if ishghandle(hwaitbar,'figure'); close(hwaitbar); end
      return;
   else
      waitbar(nstepsdone/nsteps)
   end
      
   if ~isempty(VMAPtext); displaym(VMAPtext);drawnow; end
   
end


if ~ishghandle(hwaitbar) || ~ishghandle(gcbo) || getappdata(hwaitbar,'canceled')
   if ishghandle(hwaitbar,'figure'); close(hwaitbar); end
   return;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [h,nstepsdone] = getfeature(h,hwaitbar,thislatlim,thislonlim,theme,feature,nstepsdone,nsteps)

feature = lower(feature);

% use extension of feature table filename to identify topology levels 

for i=1:length(feature)

   switch feature{i}(end-2:end)
   case 'aft'
      topolevel{i} = 'patch';
   case 'lft'
      topolevel{i} = 'line';
   case 'pft'
      topolevel{i} = 'point';
   case 'tft'
      topolevel{i} = 'text';
   otherwise
      nstepsdone = nstepsdone+1;
      waitbar(nstepsdone/nsteps)
      return
   end
   
   feature{i}(end-3:end) = [];
   
end


% get data:

% Patch

indx = strmatch('patch',topolevel);
if ~isempty(indx)
   
   for i=1:length(indx)
      
      VMAPpatch = ...
         vmap0data(h.devicename,h.library,thislatlim,thislonlim,theme,...
         'patch',feature(indx(i)));
      
      if ~ishghandle(hwaitbar) || ~ishghandle(gcbo) || getappdata(hwaitbar,'canceled')
         if ishghandle(hwaitbar,'figure'); close(hwaitbar); end
         return;
      else
         waitbar(nstepsdone/nsteps)
      end
      
      [VMAPpatch,newline,newpoint,newtext] = vmap0styles(VMAPpatch);
      
      h.extracteddata = catmlayerscell(h.extracteddata,VMAPpatch,[theme,'_',feature{indx(i)}]);
      h.extracteddata = catmlayerscell(h.extracteddata,newline  ,[theme,'_',feature{indx(i)} '_line']);
      h.extracteddata = catmlayerscell(h.extracteddata,newpoint ,[theme,'_',feature{indx(i)} '_patch']);
      h.extracteddata = catmlayerscell(h.extracteddata,newtext,[theme,'_',feature{indx(i)} '_text']);
      
      nstepsdone = nstepsdone+1/length(feature);
      
      if ~ishghandle(hwaitbar) || ~ishghandle(gcbo) || getappdata(hwaitbar,'canceled')
         if ishghandle(hwaitbar,'figure'); close(hwaitbar); end
         return;
      else
         waitbar(nstepsdone/nsteps)
      end
      
      if ~isempty(VMAPpatch); hpatch = displaym(VMAPpatch);zdatam(hpatch,-.25); end
      displaym(newline);
      displaym(newpoint);
      displaym(newtext);
   
   end
end

% Line

indx = strmatch('line',topolevel);
if ~isempty(indx)
   
   for i=1:length(indx)
      
      VMAPline = ...
         vmap0data(h.devicename,h.library,thislatlim,thislonlim,theme,...
         'line',feature(indx(i)));
      
      if ~ishghandle(hwaitbar) || ~ishghandle(gcbo) || getappdata(hwaitbar,'canceled')
         if ishghandle(hwaitbar,'figure'); close(hwaitbar); end
         return;
      else
         waitbar(nstepsdone/nsteps)
      end
      
      [VMAPline,newline,newpoint,newtext] = vmap0styles(VMAPline);
      
      h.extracteddata = catmlayerscell(h.extracteddata,VMAPline,[theme,'_',feature{indx(i)}]);
      h.extracteddata = catmlayerscell(h.extracteddata,newline  ,[theme,'_',feature{indx(i)} '_line']);
      h.extracteddata = catmlayerscell(h.extracteddata,newpoint ,[theme,'_',feature{indx(i)} '_patch']);
      h.extracteddata = catmlayerscell(h.extracteddata,newtext,[theme,'_',feature{indx(i)} '_text']);
      
      nstepsdone = nstepsdone+1/length(feature);
      
      if ~ishghandle(hwaitbar) || ~ishghandle(gcbo) || getappdata(hwaitbar,'canceled')
         if ishghandle(hwaitbar,'figure'); close(hwaitbar); end
         return;
      else
         waitbar(nstepsdone/nsteps)
      end
      
      if ~isempty(VMAPline); hline = displaym(VMAPline); end
      displaym(newline);
      displaym(newpoint);
      displaym(newtext);
      
   end
   
end

% Point

indx = strmatch('point',topolevel);
if ~isempty(indx)
   
   for i=1:length(indx)
      
      VMAPpoint = ...
         vmap0data(h.devicename,h.library,thislatlim,thislonlim,...
         theme,'point',feature(indx(i)));
      
      if ~ishghandle(hwaitbar) || ~ishghandle(gcbo) || getappdata(hwaitbar,'canceled')
         if ishghandle(hwaitbar,'figure'); close(hwaitbar); end
         return;
      else
         waitbar(nstepsdone/nsteps)
      end
      
      [VMAPpoint,newline,newpoint,newtext] = vmap0styles(VMAPpoint);
      
      h.extracteddata = catmlayerscell(h.extracteddata,VMAPpoint,[theme,'_',feature{indx(i)}]);
      h.extracteddata = catmlayerscell(h.extracteddata,newline  ,[theme,'_',feature{indx(i)} '_line']);
      h.extracteddata = catmlayerscell(h.extracteddata,newpoint ,[theme,'_',feature{indx(i)} '_patch']);
      h.extracteddata = catmlayerscell(h.extracteddata,newtext,[theme,'_',feature{indx(i)} '_text']);
      
      nstepsdone = nstepsdone+1/length(feature);
      
      if ~ishghandle(hwaitbar) || ~ishghandle(gcbo) || getappdata(hwaitbar,'canceled')
         if ishghandle(hwaitbar,'figure'); close(hwaitbar); end
         return;
      else
         waitbar(nstepsdone/nsteps)
      end
      
      if ~isempty(VMAPpoint); displaym(VMAPpoint);drawnow; end
      displaym(newline);
      displaym(newpoint);
      displaym(newtext);
      
   end
   
end

indx = strmatch('text',topolevel);
if ~isempty(indx)
   
   for i=1:length(indx)
      
      VMAPtext = ...
         vmap0data(h.devicename,h.library,thislatlim,thislonlim,...
         theme,'text',feature(indx(i)));
      
      if ~ishghandle(hwaitbar) || ~ishghandle(gcbo) || getappdata(hwaitbar,'canceled')
         if ishghandle(hwaitbar,'figure'); close(hwaitbar); end
         return;
      else
         waitbar(nstepsdone/nsteps)
      end
      
      h.extracteddata = catmlayerscell(h.extracteddata,VMAPtext,[theme,'_',feature{indx(i)}]);
      nstepsdone = nstepsdone+1/length(feature);
      
      if ~ishghandle(hwaitbar) || ~ishghandle(gcbo) || getappdata(hwaitbar,'canceled')
         if ishghandle(hwaitbar,'figure'); close(hwaitbar); end
         return;
      else
         waitbar(nstepsdone/nsteps)
      end
      
      if ~isempty(VMAPtext); displaym(VMAPtext);drawnow; end
      
   end
   
end


if ~ishghandle(hwaitbar) || ~ishghandle(gcbo) || getappdata(hwaitbar,'canceled')
   if ishghandle(hwaitbar,'figure'); close(hwaitbar); end
   return;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [h,nstepsdone] = gettheme(h,hwaitbar,thislatlim,thislonlim,theme,nstepsdone,nsteps)

% get data

[VMAPpatch,VMAPline] = ...
   vmap0data(h.devicename,h.library,thislatlim,thislonlim,theme,{'patch','line'});

nstepsdone = nstepsdone+1;
if ~ishghandle(hwaitbar) || ~ishghandle(gcbo) || getappdata(hwaitbar,'canceled')
   if ishghandle(hwaitbar,'figure'); close(hwaitbar); end
   return;
else
   waitbar(nstepsdone/nsteps)
end

[VMAPpatch,newline,newpoint,newtext] = vmap0styles(VMAPpatch);

if ~ishghandle(hwaitbar) || ~ishghandle(gcbo) || getappdata(hwaitbar,'canceled')
   if ishghandle(hwaitbar,'figure'); close(hwaitbar); end
   return;
else
   waitbar(nstepsdone/nsteps)
end

h.extracteddata = catmlayerscell(h.extracteddata,VMAPpatch,[theme 'patch']);
h.extracteddata = catmlayerscell(h.extracteddata,newline ,[theme '_line']);
h.extracteddata = catmlayerscell(h.extracteddata,newpoint,[theme '_patch']);
h.extracteddata = catmlayerscell(h.extracteddata,newtext ,[theme '_text']);

if ~isempty(VMAPpatch); hpatch = displaym(VMAPpatch);zdatam(hpatch,-.25); end
displaym(newline);
displaym(newpoint);
displaym(newtext);

nstepsdone = nstepsdone+1;
if ~ishghandle(hwaitbar) || ~ishghandle(gcbo) || getappdata(hwaitbar,'canceled')
   if ishghandle(hwaitbar,'figure'); close(hwaitbar); end
   return;
else
   waitbar(nstepsdone/nsteps)
end

[VMAPline,newline,newpoint,newtext] = vmap0styles(VMAPline);

if ~ishghandle(hwaitbar) || ~ishghandle(gcbo) || getappdata(hwaitbar,'canceled')
   if ishghandle(hwaitbar,'figure'); close(hwaitbar); end
   return;
else
   waitbar(nstepsdone/nsteps)
end

h.extracteddata = catmlayerscell(h.extracteddata,VMAPline,[theme 'line']);
h.extracteddata = catmlayerscell(h.extracteddata,newline ,[theme '_line']);
h.extracteddata = catmlayerscell(h.extracteddata,newpoint,[theme '_patch']);
h.extracteddata = catmlayerscell(h.extracteddata,newtext ,[theme '_text']);



if ~ishghandle(hwaitbar) || ~ishghandle(gcbo) || getappdata(hwaitbar,'canceled')
   if ishghandle(hwaitbar,'figure'); close(hwaitbar); end
   return;
else
   waitbar(nstepsdone/nsteps)
end

if ~isempty(VMAPline); displaym(VMAPline);drawnow; end
displaym(newline);
displaym(newpoint);
displaym(newtext);

VMAPpoint = ...
   vmap0data(h.devicename,h.library,thislatlim,thislonlim,theme,'point');

nstepsdone = nstepsdone+1;
if ~ishghandle(hwaitbar) || ~ishghandle(gcbo) || getappdata(hwaitbar,'canceled')
   if ishghandle(hwaitbar,'figure'); close(hwaitbar); end
   return;
else
   waitbar(nstepsdone/nsteps)
end
   
[VMAPpoint,newline,newpoint,newtext] = vmap0styles(VMAPpoint);

h.extracteddata = catmlayerscell(h.extracteddata,VMAPpoint,[theme 'point']);
h.extracteddata = catmlayerscell(h.extracteddata,newline ,[theme '_line']);
h.extracteddata = catmlayerscell(h.extracteddata,newpoint,[theme '_patch']);
h.extracteddata = catmlayerscell(h.extracteddata,newtext ,[theme '_text']);

if ~ishghandle(hwaitbar) || ~ishghandle(gcbo) || getappdata(hwaitbar,'canceled')
   if ishghandle(hwaitbar,'figure'); close(hwaitbar); end
   return;
else
   waitbar(nstepsdone/nsteps)
end

if ~isempty(VMAPpoint); displaym(VMAPpoint);drawnow; end
displaym(newline);
displaym(newpoint);
displaym(newtext);

VMAPtext = ...
   vmap0data(h.devicename,h.library,thislatlim,thislonlim,theme,'text');

nstepsdone = nstepsdone+1;
if ~ishghandle(hwaitbar) || ~ishghandle(gcbo) || getappdata(hwaitbar,'canceled')
   if ishghandle(hwaitbar,'figure'); close(hwaitbar); end
   return;
else
   waitbar(nstepsdone/nsteps)
end

h.extracteddata = catmlayerscell(h.extracteddata,VMAPtext,[theme 'text']);
if ~isempty(VMAPtext); displaym(VMAPtext);drawnow; end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [theme,feature,propval,nsteps]=assembleRequests(indx,n,VMAPnames,objhndl,objparent,h)

theme = [];
feature = [];
propval = [];
nsteps = 0;

k=0;
for i=1:length(n)
   
   switch n(i)
   case 0
      
      k=k+1;
      theme{k} = VMAPnames{i}.coverage_name;
      feature{k} = [];
      propval{k} = [];
      nsteps = nsteps+4;
      
   case 1
      
      thisfeaturename = lower(VMAPnames{i}.name);
      thisthemename = h.VMAPnames{ ...
            find( objhndl == objparent(indx(i)) ) ...
         }.coverage_name;
      indxtheme = strmatch(thisthemename,theme);
      for j=length(indxtheme):-1:1
         if ~isempty(propval{indxtheme(j)}); indxtheme(j) = []; end
      end
      
      if isempty(indxtheme) 
         k=k+1;
         nsteps = nsteps + 1;
      	  theme{k} = thisthemename;
         feature{k} = {thisfeaturename};
         propval{k} = [];
      else
         fcell = feature{indxtheme};
         feature{indxtheme} = {fcell{:},thisfeaturename};
      end
      
   case 2
      
      thisthemename = h.VMAPnames{...
            find( objhndl == ...
               objparent( ...
                  find( objhndl == objparent(indx(i)) ) ...
               ) ...   
               ) ...
         }.coverage_name;
      
      thisfeaturename = ...
         lower( ...
            h.VMAPnames{...
                  find( objhndl == objparent(indx(i)) ) ...
               }.name ...
            );
         
         indxtheme = strmatch(thisthemename,theme);
         indxfeature = [];
         for j=1:length(indxtheme)
            indxfeature = strmatch(thisfeaturename,feature{indxtheme(j)});
            if ~isempty(indxfeature); 
               indxtheme = indxtheme(j);
               break; 
            end
         end
         
         if isempty(indxfeature)
            k=k+1;
            nsteps = nsteps + 1;
      	     theme{k} = thisthemename;
            feature{k} = {thisfeaturename};
            propval{k} = {VMAPnames{i}.attribute,VMAPnames{i}.value};
         else
            propval{k} = {propval{k}{:},VMAPnames{i}.attribute,VMAPnames{i}.value};
         end
         
      
   end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function mcell = catmlayerscell(mcell,struc,name)

% CATMLAYERSCELL concatenates a structure to a MLAYERS-style cell array.
% This cell array is N by 2 in size, with the first column containing 
% display structures, while the second contains associated 
% variable names. If the structure to be added matches an existing name,
% the structure is merged with that existing structure. Otherwise, a 
% new entry is added to the cell array.

if isempty(mcell) && isempty(struc)
   return
elseif isempty(mcell) && ~isempty(struc)
   mcell = {struc,name};
   return
elseif ~isempty(mcell) && isempty(struc)
   return
end

indx = strmatch(name,mcell(:,2),'exact');

if isempty(indx) && ~isempty(struc)
   mcell{end+1,1} = struc;
   mcell{end,2} = name;
elseif ~isempty(struc)
   mcell{indx(1),1} = catstructures(mcell{indx(1),1},struc);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function s1 = catstructures(s1,s2)

start = length(s1);
if start ==0
    s1 = s2;
elseif isempty(s2)
    s1 = s1;
else
    names = union(fieldnames(s1),fieldnames(s2));

% enforce same fieldnames for each structure    

    l1 = length(s1);
    l2 = length(s2);
    for i=1:length(names)
        s1 = setfield(s1,{l1+1},names{i},[]);
        s2 = setfield(s2,{l2+1},names{i},[]);
    end
    s1(l1+1) = [];
    s2(l2+1) = [];
    
% convert structures to cells
    [ignored,indx1] = sort(fieldnames(s1));
    cell1 = struct2cell(s1);
    cell1 = cell1(indx1,:,:);
    
    [ignored,indx2] = sort(fieldnames(s2));
    cell2 = struct2cell(s2);
    cell2 = cell2(indx2,:,:);

% concatenate cells

    cell1(:,:,l1 + (1:l2)) = deal(cell2(:,:,1:l2));
    
    s1 = cell2struct(cell1,names,1);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function vmap0uihelp(~,~)
% Displays help text for VMAP0UI

%check if a get is in progress
if ~isempty( findall(0,'type','figure','tag','TMWWaitbar','UserData',get(gcbo,'Parent')) );
   return
end

helpstr = { ...
      [],...
      [],'OVERVIEW',[], ...
      [ ...
         'The VMAP0UI screen lets you read data from the Vector Map Level 0 (VMAP0). ' ...
         'The VMAP0 is the most detailed map database available to the public. ' ...
      ],[],[ ...
         'You use the list to select the type of data and the map to select the region of interest. ' ...
         'When you click the GET button, data is extracted and displayed on the map. ' ...
         'Use the SAVE button to save the data in a MAT-file or to the base workspace for later display. ' ...
         'The CLOSE button closes the window.' ...
         '' ...
         '' ...
         '' ...
      ],[],[],'THE MAP',[],[ ...
         'The map controls the geographic extent of the data to be extracted. ' ...
         'VMAP0UI extracts data for areas currently visible on the map. ' ...
         'Use the mouse to zoom in or out to the area of interest. ' ...
         'Type HELP ZOOM for more on zooming. ' ...
      ],[],[ ...
         'The VMAP0 divides the world up into tiles of about 5 by 5 degrees. ' ...
         'When extracting, data is returned for all visible tiles, including those parts of the tile that are outside of the current view. ' ...
         'The map shows the VMAP0 tiles in light yellow with light gray edges. ' ...
         'The data density is high, so extracting data for a large number of tiles may take much time and memory.' ...
         'A count of the number of visible tiles is above the map. ' ...
      ],[],[],'THE LIST',[],[ ...
         'The list controls the type of data to be extracted. ' ...
         'The tree structure of the list reflects the structure of the VMAP0 database. ' ... 
         'Upon starting VMAP0UI, the list shows the major categories of VMAP data, called themes.' ...
         'Themes are subdivided into features, which consist of data of common graphic type (patch, line, point or text) or cultural type (airport, roads, railroads). ' ...
         'Double-click on a theme to see the associated features. ' ...
         'Features may have properties and values, for example a railroad tracks property, with values single or multiple. ' ...
         'Double-click on a feature to see the associated properties and values. ' ...
         'Double clicking on an open theme or feature will close it. ' ...
         'When a theme is selected, VMAP0UI will get all of the associated features. ' ...
         'When a feature is selected, VMAP0UI will get all of that feature data. ' ...
         'When properties and values are selected, VMAP0UI will get the data for any of the properties and values match (i.e., the union operation). ' ... 
      ],[],[],'GET BUTTON',[],[ ...
         'The GET button reads the currently selected VMAP0 data and displays it on the map. ' ...
         'Use the CANCEL button on the progress bar to interrupt the process. ' ...
         'For a quicker response, press the standard interrupt key combination for your platform. ' ...
      ],[],[],'CLEAR BUTTON',[],[ ...
         'The CLEAR button removes any previously read VMAP0 data from the map. ' ...
      ],[],[],'SAVE BUTTON',[],[ ...
         'The SAVE button saves the currently displayed VMAP0 data to a Mat-file or the base workspace. ' ...
         'If you choose to save to a file, you will be prompted for a file name and location. ' ...
         'If you choose to save to the base workspace, you will be notified of the variable names that will be overwritten. ' ...
         'The results are stored as display structures with variables names based on theme and feature names. ' ...
         'Use LOAD and DISPLAYM to redisplay the data from a file on a map axes. ' ...
         'You can also use the MLAYERS GUI to read and display the data from a file. ' ...
         'To display the data in the base workspace, use DISPLAYM. ' ...
         'To display all of the display structures, use ROOTLAYR; DISPLAYM(ans). ' ...
         'To display all of the display structures using the MLAYERS GUI, type ROOTLAYR; MLAYERS(ans). ' ...
      ],[],[],'CLOSE BUTTON',[],[ ...
         'The CLOSE button closes the VMAP0 panel ' ...
      ]} ...
   ;


h.fig = figure('Visible','on','Units','character','Position',[5 5 80 40],'Resize','off','Menubar','none','Name','VMAP0UI Help','NumberTitle','off');
h.list = uicontrol('Style','list','units','normalized','position',[.1 .15 .7 .8],'BackgroundColor','w','SelectionHighlight','off');
helpstr = textwrap(h.list,helpstr);

% make room for scrollbar
set(h.list,'String',helpstr,'Position',[.1 .15 .8 .8],'ListboxTop',2)
h.button = uicontrol('Style','push','Units','Character','Position',[35 2 10 2],'String','Close','Callback','close(gcbf)');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function vmap0uilist(~,~)
% Manages insertions and deletions to the list of themes, features and
% property-value pairs

h = get(gcbf,'Userdata');

switch get(gcbf,'SelectionType')
case 'open'
    
    value = get(gcbo,'Value');
    strs = get(gcbo,'string');
    
    % do operation on top level objects only
    ph = h.objhndl(value);
    gens = gennum(ph,h.objhndl,h.objparent);
    value = value(gens==min(gens));
    
    selected = value;
    
    for i=length(value):-1:1
        
       switch gennum( h.objhndl(value(i)) ,h.objhndl,h.objparent)+1
       
       case 1 % 'theme'

            switch h.expandable(value(i))
            case 0                            % remove subclasses

                strs(value(i),2) = '+';
                
                h.expandable(value(i)) = 1;
                
                chil = allekinder(h.objhndl(value(i)),h.objhndl,h.objparent);
                rowsToDelete = find(ismember(h.objhndl,chil));

                strs(rowsToDelete,:) = [];
                selected = value(i);
                
                h.objhndl(rowsToDelete) = [];
                h.objparent(rowsToDelete) = [];
                h.VMAPnames(rowsToDelete) = [];
                h.expandable(rowsToDelete) = [];

            case 1                            % insert subclasses
               
                h.expandable(value(i)) = 0;
                
                suffixstr = '     +   ';
                names = str2mat(h.FTfilenames{h.objhndl(value(i))}.description);
                newstrs = [repmat(suffixstr,size(names,1),1) names];
                strs = str2mat( ...
                    strs(1:value(i),:) , ...
                    newstrs , ...
                    strs(value(i)+1:end,:) ...
                    );
                    
                pseudohandle = h.objhndl(value(i));
                hvec = h.objhndl;
                pvec = h.objparent;
                [hvec,pvec] = insertkinder(pseudohandle,size(names,1),hvec,pvec);
                
                h.objhndl = hvec ;
                h.objparent = pvec;
                
                struct2insert = h.FTfilenames{h.objhndl(value(i))};
                for ii=1:length(struct2insert)
                   cells2insert{ii} = struct2insert(ii);
                end   
                h.VMAPnames = { ...
                      h.VMAPnames{1:value(i)},...
                      cells2insert{:},...
                      h.VMAPnames{value(i)+1:end} ...
                   };
                
                h.expandable = [ ...
                      h.expandable(1:value(i)),...
                      ones(1,length(cells2insert)),...
                      h.expandable(value(i)+1:end) ...
                   ];
                
           end % switch on whether expandable
       
           set(gcbo,'String',strs,'Value',selected)
           
       case 2 % 'Feature'

            switch h.expandable(value(i)) % strs(value(i),2)
            case 0                            % remove subclasses

                strs(value(i),6) = '+';
                
                h.expandable(value(i)) = 1;
                
                chil = allekinder(h.objhndl(value(i)),h.objhndl,h.objparent);
                rowsToDelete = find(ismember(h.objhndl,chil));

                strs(rowsToDelete,:) = [];
                selected = value(i);
                
                h.objhndl(rowsToDelete) = [];
                h.objparent(rowsToDelete) = [];
                h.VMAPnames(rowsToDelete) = [];
                h.expandable(rowsToDelete) = [];
                
             case 1                            % add subclasses
                
                strs(value(i),6) = '-';
                
                h.expandable(value(i)) = 0;
                
                themenum = h.objparent(value(i));
                featnum = sum(h.objparent(h.objparent(value(i)):value(i))== h.objparent(value(i)));
                
                FT = h.FTfilenames{themenum};
                
                featuretablename = lower(FT(featnum).name);
                FTattributes = {FT(featnum).fields.name};
                ivdt = h.ivdt{themenum};
                cvdt = h.cvdt{themenum};
                                    
                indx = [];
                if ~isempty(ivdt); 
                    indx = strmatch(featuretablename,{ivdt.table});
                end
         
                if ~isempty(indx)
                   ivdtnames = str2mat(ivdt(indx).description);

                    
                    iattributeDescription = {};
                    for j=1:length(ivdt(indx))
                        thisvdt = ivdt(indx(j));
                        [junk,FTindx] = intersect(FTattributes,{thisvdt.attribute});
                        iattributeDescription{j} = FT(featnum).fields(FTindx).description;
                    end
                    iattributeDescription = str2mat(iattributeDescription{:});
                        
                    
                else
                   ivdtnames = '';
                   iattributeDescription = '';
                end
                iindx = indx; 
                
                
                indx = [];
                if ~isempty(cvdt); 
                    indx = strmatch(featuretablename,{cvdt.table});
                end
                
                if ~isempty(indx)
                    cvdtnames = str2mat(cvdt(indx).description);
                                        
                    cattributeDescription = {};
                    for j=1:length(cvdt(indx))
                        thisvdt = cvdt(indx(j));
                        [junk,FTindx] = intersect(FTattributes,{thisvdt.attribute});
                        cattributeDescription{j} = FT(featnum).fields(FTindx).description;
                    end
                    cattributeDescription = str2mat(cattributeDescription{:});
                        
                else
                   cvdtnames = '';
                   cattributeDescription = '';
                end
                cindx = indx;
                
                names = strvcat(ivdtnames,cvdtnames);
                attdescriptions = strvcat(iattributeDescription,cattributeDescription);
                              
                suffixstr = '             ';
                
                nrows = size(names,1);
                newstrs = [repmat(suffixstr,nrows,1) attdescriptions repmat(': ',nrows,1) names];
                
                if nrows > 0
                   strs = str2mat( ...
                      strs(1:value(i),:) , ...
                      newstrs , ...
                      strs(value(i)+1:end,:) ...
                      );
                else
                   % nothing to expand, so remove indication that there might be more 
                   strs(value(i),2) = ' ';                 
                end
                
                pseudohandle = h.objhndl(value(i));
                hvec = h.objhndl;
                pvec = h.objparent;
                    
                [hvec,pvec] = insertkinder(pseudohandle,size(names,1),hvec,pvec);
                
                h.objhndl = hvec;
                h.objparent = pvec;
                
                cells2insert = [];
                for ii=1:length(iindx)
                   cells2insert{ii} = ivdt(iindx(ii));
                end   
                for ii=1:length(cindx)
                   cells2insert{ii+length(iindx)} = cvdt(cindx(ii));
                end   
                
                if ~isempty(cells2insert)
                   h.VMAPnames = { ...
                         h.VMAPnames{1:value(i)},...
                         cells2insert{:},...
                         h.VMAPnames{value(i)+1:end} ...
                      };
                end
                
                h.expandable = [ ...
                      h.expandable(1:value(i)),...
                      zeros(1,length(cells2insert)),...
                      h.expandable(value(i)+1:end) ...
                   ];

                

            end % switch on +/- (whether expanded)
        
           set(gcbo,'String',strs,'Value',selected)

        end % switch on level
    
    end % for


end % switch selection type

set(gcbf,'UserData',h)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function hc = allekinder(h,hvec,pvec)
% ALLEKINDER descendents of an object based on pseudohandles
%
% hc = ALLEKINDER(h,hvec,pvec) finds the descendents of an object with 
% pseudohandle h given the vector of pseudohandles hvec and the associated 
% vector of parent pseudohandles pvec.
%
% Pseudohandles are unique numbers associated with a set of heirarchical 
% objects. They are tracked by means of the object and parent vectors hvec 
% and pvec. For example, three objects, the second of which has two
% subobjects, might be described by the following vectors:
%
%      hvec = [ 1 2 4 5 3] and pvec = [ 0 0 2 2 0]
% 
% The pseudohandles of the children of the object with the pseudohandle 2
% are then [4 5]. 

hc = [];

while 1
    
    thesekids = hvec(find(ismember(pvec,h)));
    if isempty(thesekids); break; end;
    hc = [hc thesekids];
    h = thesekids;

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function n = gennum(h,hvec,pvec)
% GENNUM generation number of an object based on pseudohandles
%
% n = GENNUM(h,hvec,pvec) finds the generation of an object with 
% pseudohandle h given the vector of pseudohandles hvec and the associated 
% vector of parent pseudohandles pvec. Objects that have no parents have
% a generation number of 0, while the grandchildren of such objects have a
% generation number of 2
%
% Pseudohandles are unique numbers associated with a set of heirarchical 
% objects. They are tracked by means of the object and parent vectors hvec 
% and pvec. For example, three objects, the second of which has two subobjects, 
% might be described by the following vectors:
%
%      hvec = [ 1 2 4 5 3] and pvec = [ 0 0 2 2 0]
% 
% The pseudohandles of the children of the object with the pseudohandle 2 are 
% then [4 5]. 

n = -ones(size(h));

for i=1:length(h)
    while 1
        
        indx = find(ismember(hvec,h(i)));
        theseparents = pvec(indx);
        if isempty(theseparents); break; end
        h(i) = unique(theseparents);
        n(i) = n(i)+1;
    
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [hvec,pvec] = insertkinder(h,nchildren,hvec,pvec)

if length(h) > 1
    error('map:vmap0ui:nonscalarHandle', 'One object at a time, please.')
end

indx = find(h==hvec);

hvec = hvec(:)';
pvec = pvec(:)';


hchildren = max(hvec) + (1:nchildren);
hparents  = h * ones(size(hchildren));

pvec = ...
    [ ...
    pvec(1:indx)  ...
    hparents  ...
    pvec(indx+1:end) ...
    ];
 
hvec = [ ...
    hvec(1:indx)  ...
    hchildren  ...
    hvec(indx+1:end) ...
    ];
