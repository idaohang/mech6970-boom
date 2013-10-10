function [map,refvec,UHL,DSI,ACC] = dted(varargin)
%DTED Read U.S. Dept. of Defense Digital Terrain Elevation Data (DTED)
%
%   [Z, REFVEC] = DTED returns all of the elevation data in a DTED file as
%   a regular data grid with elevations in meters.  The file is selected
%   interactively.  This function reads the DTED elevation files, which
%   generally have filenames ending in ".dtN", where N is 0,1,2,3,...
%   REFVEC is the associated referencing vector.
%
%   [Z, REFVEC] = DTED(FILENAME) returns all of the elevation data in the
%   specified DTED file.  The file must be found on the MATLAB path. If not
%   found, the file may be selected interactively.
%
%   [Z, REFVEC] = DTED(FILENAME, SAMPLEFACTOR) subsamples data from the 
%   specified DTED file.  SAMPLEFACTOR is a scalar integer.  When
%   SAMPLEFACTOR is 1 (the default), DTED reads the data at its full
%   resolution.  When SAMPLEFACTOR is an integer n greater than one, every
%   n-th point is read.
% 
%   [Z, REFVEC] = DTED(FILENAME, SAMPLEFACTOR, LATLIM, LONLIM) reads the
%   data  for the part of the DTED file within the latitude and longitude
%   limits.  The limits must be two-element vectors in units of degrees. 
%
%   [Z, REFVEC] = DTED(DIRNAME, SAMPLEFACTOR, LATLIM, LONLIM) reads and 
%   concatenates data from multiple files within a DTED CD-ROM or directory 
%   structure.  The dirname input is a string with the name of a directory 
%   containing the DTED directory.  Within the DTED directory are
%   subdirectories for each degree of longitude, each of which contain
%   files for each degree of latitude.  For DTED CD-ROMs, dirname is the
%   device name of the CD-ROM drive.  LATLIM may not span 50 degrees
%   North or South latitude.
%
%   [Z, REFVEC, UHL, DSI, ACC] = DTED(...) returns structures containing
%   the DTED User Header Label (UHL), Data Set Identification (DSI) and
%   ACCuracy metadata records.
%
%   Latitude-Dependent Sampling
%   ---------------------------
%   In DTED files north of 50 degrees North and south of 50 degrees
%   South, where the meridians have converged significantly relative to
%   the equator, the longitude sampling interval is increased to twice
%   the latitude sampling interval, from 30" to 60" in the case of Level
%   0, for example.  In this case, in order to retain square output
%   cells, this function changes the latitude sampling to match the
%   longitude sampling. For example, it will return a 121-by-121
%   elevation grid for a DTED file covering 49 to 50 degrees north, but
%   a 61-by-61 grid for a file covering 50 to 51 degrees north.
%
%   If a directory name is supplied instead of a file name and LATLIM
%   spans either 50 degrees North or 50 degrees South, an error results.
%
%     LATLIM = [20 60];   <-- error    LATLIM = [-55 -45];   <-- error
%     LATLIM = [50 60];   <-- OK       LATLIM = [-55 -50];   <-- OK
%     LATLIM = [20 50];   <-- OK       LATLIM = [-50 -45];   <-- OK
%
%   Going north the longitude sampling interval increases further at the
%   latitudes of 70 degrees N, 75 N, and 80 N.  Likewise, going south,
%   there is a increase at 70 S, 75 S, and 80 S.  LATLIM must not span
%   any of these latitudes either.
%
%   Null Data Values
%   ----------------
%   Some DTED Level 1 and higher data tiles contain null data tiles, coded
%   with value -32767.  When encountered, these null data values are
%   converted to NaN.
%
%   Non-Conforming Data Encoding
%   ----------------------------
%   DTED files from some sources may depart from the specification by using
%   twos-complement encoding for binary elevation files instead of
%   "sign-bit" encoding.  This difference affects the decoding of negative
%   values, and incorrect decoding usually leads to nonsensical elevations.
%   Thus, if the DTED function determines that all the (non-null) negative
%   values in a file would otherwise be less than -12,000 meters, it issues
%   a warning and assumes twos-complement encoding.
%
%   Data Sources and Information
%   ----------------------------
%   DTED files contain digital elevation maps covering 1-by-1-degree 
%   quadrangles at horizontal resolutions ranging from about 1 kilometer to
%   1 meter.  For details on locating DTED for download over the Internet,
%   see the following documentation at the MathWorks web site: 
% 
%   <a href="matlab: 
%   web('http://www.mathworks.com/support/tech-notes/2100/2101.html#dted') 
%   ">http://www.mathworks.com/support/tech-notes/2100/2101.html</a>
%
%   See also DTEDS.

% Copyright 1996-2008 The MathWorks, Inc.
% $Revision: 1.1.6.16 $  $Date: 2009/01/16 11:04:51 $

if nargin < 1; 
    [map,refvec,UHL,DSI,ACC] = dtedf; 
    return
end

name = varargin{1};
if exist(name,'dir') == 7
    if nargin < 4
        error(['map:' mfilename ':mapformatsError'], ...
            'Latlim and lonlim required for directory calling form')
    end
    [map,refvec,UHL,DSI,ACC] = dtedc(varargin{:});
else
    [map,refvec,UHL,DSI,ACC] = dtedf(varargin{:});
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [map,refvec,UHL,DSI,ACC] = dtedc(rd,samplefactor,latlim,lonlim)

% Concatenate adjacent DTED tiles, accounting for the fact that the
% edges of adjacent tiles contain redundant data.

% add file separator to root directory if necessary
if ~strcmp(rd(end),filesep)
    rd(end+1) = filesep;
end

% error checking for input arguments
error(nargchk(1,4,nargin,'struct'))

if nargin < 2
    samplefactor = 1;
end

if nargin < 3
    latlim = [];
end

if nargin < 4
    lonlim = [];
end

% If request just touches edge of the next tile, don't read it
if mod(latlim(2),1) == 0
    latlim(2) = latlim(2)-epsm('deg');
end

if mod(lonlim(2),1) == 0;
    lonlim(2) = lonlim(2)-epsm('deg');
end

% round the limits since DTED is read in 1 deg x 1 deg square tiles
latmin = floor(latlim(1));
latmax = floor(latlim(2));
lonmin = floor(lonlim(1));
lonmax = floor(lonlim(2));

% LATLIM must not span +/- 50 degrees
if (latmin < -50 && latmax > -50) || (latmin < 50 && latmax > 50)
    eid = sprintf('%s:%s:latlimSpans50',getcomp,mfilename); 
    error(eid,'LATLIM spans 50 degrees North or South')
end

% define columns and rows for tiles to be read
uniquelons = lonmin:lonmax;
uniquelats = latmin:latmax;
dtedfile = cell(numel(uniquelats), numel(uniquelons));
levels   = cell(numel(uniquelats), numel(uniquelons));

% redefine uniquelons if lonlim extends across the International Dateline
if lonmin > lonmax
    indx1 = lonmin:179;
    indx2 = -180:lonmax;
    uniquelons = [indx1 indx2];
end

[latdir,londir] = dteddirs(uniquelons,uniquelats,'');

% check to see if the files exist
for k = 1:length(uniquelats)
    for j = 1:length(uniquelons)
        for n = 0:3;
            filename = [rd londir{j} latdir{k} 'dt' num2str(n)];
            if exist(filename,'file') == 2
                dtedfile{k,j} = filename;
                levels{k,j} = n;
                break
            end
            filename(end) = '*';
            dtedfile{k,j} = filename;
        end
    end
end

% trim off requests for missing files around edges
changed = 1;
while changed
   changed = 0;
   if ~isempty(levels) && isempty([ levels{:,1} ])
      dtedfile(:,1) = [];
      levels(:,1) = [];
      changed = 1;
   end
   if ~isempty(levels) && isempty([ levels{:,end} ])
      dtedfile(:,end) = [];
      levels(:,end) = [];
      changed = 1;
   end
   if ~isempty(levels) && isempty([ levels{1,:} ])
      dtedfile(1,:) = [];
      levels(1,:) = [];
      changed = 1;
   end
   if ~isempty(levels) && isempty([ levels{end,:} ])
      dtedfile(end,:) = [];
      levels(end,:) = [];
      changed = 1;
   end
end

% Stop if missing files
if isempty(dtedfile)
    eid = sprintf('%s:%s:noData',getcomp,mfilename);
    error(eid,'No data for requested area.')
end

% break out if only 1 tile is required
if numel(dtedfile) == 1
    [map,refvec,UHL,DSI,ACC] = ...
        dted(dtedfile{1,1},samplefactor,latlim,lonlim);
   return
end

level = unique([levels{:}]);
if length(level)>1
    error(['map:' mfilename ':mapformatsError'], ...
        'Inconsistent Levels between files. Check path or CD.')
end

nrowMat = NaN(size(dtedfile));  
ncolMat = nrowMat;
% read all files to compute number of rows and number of columns in each tile
for k = 1:size(dtedfile,1) 
    for j = 1:size(dtedfile,2)
        if exist(dtedfile{k,j}, 'file') == 2
            tmap = dted(dtedfile{k,j},samplefactor,latlim,lonlim);  
            nrowMat(k,j) = size(tmap,1);
            ncolMat(k,j) = size(tmap,2);
        end
    end
end

% replace nans with the values required for correct concatenation
nrows = max(nrowMat,[],2);  nrows(isnan(nrows)) = max(nrows);
ncols = max(ncolMat,[],1);  ncols(isnan(ncols)) = max(ncols);

map = cell(size(dtedfile,1),1);
refvec = cell(size(dtedfile,1),1);

% read the first file (bottom left hand corner of map)
if exist(dtedfile{1,1}, 'file') == 2
    [map{1},refvec{1},UHL(1,1),DSI(1,1),ACC(1,1)] = ...
        dted(dtedfile{1,1},samplefactor,latlim,lonlim);
else
    % If the first file does not exist, determine what the refvec would
    % be if it did exist.  Also assign metadata structures just to get the
    % structure field names.
    [map{1},refvec{1},UHL(1,1),DSI(1,1),ACC(1,1)] =...
        readFirstNonEmpty(dtedfile,samplefactor,latlim,lonlim);
    map{1} = NaN(nrows(1), ncols(1));
end

% Create structures with fields that contain no data. We'll use it if we're
% missing a data file.
UHLempty = UHL(1,1);
fdnames = fieldnames(UHLempty);
for i = 1:length(fdnames)
    UHLempty.(fdnames{i}) = '';
end
DSIempty = DSI(1,1);
fdnames = fieldnames(DSIempty);
for i = 1:length(fdnames)
    DSIempty.(fdnames{i}) = '';
end
ACCempty = ACC(1,1);
fdnames = fieldnames(ACCempty);
for i = 1:length(fdnames)
    ACCempty.(fdnames{i}) = '';
end

if size(dtedfile,1) > 1
    % read remaining files in 1-st column (same longitude, increasing latitudes)
    for k = 2:size(dtedfile,1)
        if exist(dtedfile{k,1}, 'file') == 2
            [map{k},refvec{k},UHL(k,1),DSI(k,1),ACC(k,1)] ...
                = dted(dtedfile{k,1},samplefactor,latlim,lonlim);
        else
            map{k} = NaN(nrows(k), ncols(1));
            refvec{k} = refvec{k-1} ...
                + [0 (nrows(k)-1) 0]/refvec{k-1}(1);
            UHL(k,1) = UHLempty;
            DSI(k,1) = DSIempty;
            ACC(k,1) = ACCempty;
        end
        % Strip the first row off each tile, because it's
        % already contained in the preceding tile.
        map{k}(1,:) = [];
    end
    frefvec = refvec{k};
else
    frefvec = refvec{1};
end

fmap = cell(size(dtedfile,2),1);

% concatenate tiles in the first column
fmap{1} = cat(1,map{:});

% read remaining files for remaining columns and rows
% do one column on each pass through the outer loop.
for j = 2:size(dtedfile,2)
    
    % Clear out each cell so we can re-use the map and refvec arrays
    [map{:}] = deal([]);
    [refvec{:}] = deal([]);
    
    if exist(dtedfile{1,j}, 'file') ~= 0
        % read file corresponding to the 1-st (bottom) tile in the j-th column
        [map{1},refvec{1},UHL(1,j),DSI(1,j),ACC(1,j)] ...
            = dted(dtedfile{1,j},samplefactor,latlim,lonlim);
    else
        map{1} = NaN(nrows(1), ncols(j));
        UHL(1,j) = UHLempty;
        DSI(1,j) = DSIempty;
        ACC(1,j) = ACCempty;
    end
    % Strip the first column off each tile, because it's
    % already contained preceding column of tiles.
    map{1}(:,1) = [];
    % read the remaining files in the j-th column
    for k = 2:size(dtedfile,1)
        if exist(dtedfile{k,j}, 'file') == 2
            [map{k},refvec{k},UHL(k,j),DSI(k,j),ACC(k,j)] ...
                = dted(dtedfile{k,j},samplefactor,latlim,lonlim);
        else
            map{k} = NaN(nrows(k), ncols(j));
            UHL(k,j) = UHLempty;
            DSI(k,j) = DSIempty;
            ACC(k,j) = ACCempty;
        end
        map{k}(1,:) = [];  % Strip off the first row
        map{k}(:,1) = [];  % Strip off the first column
    end
    % concatenate the tiles in the j-th column
    fmap{j} = cat(1,map{:});
end

% concatenate the columns
map = cat(2,fmap{:});
refvec = frefvec;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [map,refvec,UHL,DSI,ACC] = dtedf(filename,scalefactor,latlim,lonlim)

if nargin==0;
   filename = [];
   scalefactor = 1;
   latlim = [];
   lonlim = [];
elseif nargin==1
   scalefactor = 1;
   latlim = [];
   lonlim = [];
elseif nargin==2
   latlim = [];
   lonlim = [];
elseif nargin==3
   lonlim = [];
elseif nargin ~=4
   error(['map:' mfilename ':mapformatsError'], 'Incorrect number of arguments')
end

% ensure row vectors

latlim = latlim(:)';
lonlim = npi2pi(lonlim(:)'); % No effort made (yet) to work across the dateline
if ~isempty(lonlim) && (lonlim(2) < lonlim(1))
    lonlim(2) = lonlim(2) + 360;
end

% check input arguments

if ~isempty(filename) && ~ischar(filename)
    error(['map:' mfilename ':mapformatsError'], 'Filename must be a string');
end

if length(scalefactor) > 1
    error(['map:' mfilename ':mapformatsError'], 'Scalefactor must be a scalar');
end

if ~isempty(latlim) && ~isequal([1 2],size(latlim));
    error(['map:' mfilename ':mapformatsError'], 'latlim must be a two element vector in units of degrees');
end 

if ~isempty(lonlim) && ~isequal([1 2],size(lonlim));
    error(['map:' mfilename ':mapformatsError'], 'lonlim must be a two element vector in units of degrees');
end 

% Open the file to read the header information

fid = -1;
if ~ isempty(filename); 
   fid = fopen(filename,'rb','ieee-be');
end

if fid==-1
   [filename, path] = uigetfile('*.*', 'Please select the DTED file');
   if filename == 0
       map = [];
       refvec = [];
       UHL = [];
       DSI = [];
       ACC = [];
       return;
   end
   filename = [path filename];
   fid = fopen(filename,'rb','ieee-be');
end

% Define the header format, and read the header records

UHLfield = UHLdescription;
DSIfield = DSIdescription;
ACCfield = ACCdescription;

UHL = readfields(filename,UHLfield,1,'ieee-be',fid);
DSI = readfields(filename,DSIfield,1,'ieee-be',fid);
ACC = readfields(filename,ACCfield,1,'ieee-be',fid);

fclose(fid);

% Data records information
%
% True longitude = longitude count x data interval + origin (Offset from the SW corner longitude)
% 
% True latitude = latitude count x data interval + origin (Offset from the SW corner latitude)
%
% 1x1 degree tiles, including all edges, edges duplicated across files.


% special correction for incorrectly formatted data files for region just north
% of equator (LSJ 01182003)
if ~isempty(strfind(filename,'n00.dt0'))
    if isequal(DSI.LatitudeofNWcorner,'010000S')
        % fix data (latitude of northern points should be at 1 deg North
        % not 1 deg South)
        DSI.LatitudeofNWcorner = '010000N';
        DSI.LatitudeofNEcorner = DSI.LatitudeofNWcorner;
    end
end

% Detect and correct errors in DTED data files for cells bordered on the
% west by the prime meridian.  The longitudes found in the file may be
% inconsistent with the directory in which the data file is located.  In
% this case, assume that the directory is correct.
[pathstr, name, ext] = fileparts(filename);
[pathstr, nameOfEnclosingDirectory] = fileparts(pathstr);
if strcmp(nameOfEnclosingDirectory,'e000')
    if 'W' == DSI.LongitudeofNEcorner(end)
        % Replace 'W' with 'E' in longitude corner/origin strings.
        wid = sprintf('%s:%s:correctingLongitudes', getcomp, mfilename);
        warning(wid,...
            'Correcting ''W'' to ''E'' in DSI & UHL longitude strings for %s.',...
            [fullfile(nameOfEnclosingDirectory,name) ext]);
        DSI.Longitudeoforigin(end) = 'E';
        UHL.Longitudeoforigin(end) = 'E';
        DSI.LongitudeofSWcorner(end) = 'E';
        DSI.LongitudeofNWcorner(end) = 'E';
        DSI.LongitudeofNEcorner(end) = 'E';
        DSI.LongitudeofSEcorner(end) = 'E';
    end
end

maplatlim = [ dmsstr2deg(DSI.LatitudeofSWcorner)  dmsstr2deg(DSI.LatitudeofNWcorner) ];
maplonlim = [ dmsstr2deg(DSI.LongitudeofSWcorner) dmsstr2deg(DSI.LongitudeofSEcorner) ];
dlat = secstr2deg(DSI.Latitudeinterval);
dlon = secstr2deg(DSI.Longitudeinterval);
ncols = round(diff(maplatlim/dlat)) + 1;
nrows = round(diff(maplonlim/dlon)) + 1;

lato = dmsstr2deg(DSI.Latitudeoforigin);
lono = dmsstr2deg(DSI.Longitudeoforigin);

skipfactor = 1;
[dlat0,dlon0] = deal(dlat,dlon);
if dlat ~= dlon
   wid = sprintf('%s:%s:ignoringLatitudeSpacing', getcomp, mfilename);
   warning(wid,'%s','Latitude and longitude spacing differ; Using coarser grid')
   skipfactor = dlon/dlat;
   dlat = max([dlat dlon]);
   dlon = dlat;
end

%  Check to see if latlim and lonlim within map limits

if isempty(latlim); latlim = maplatlim; end
if isempty(lonlim); lonlim = maplonlim; end


errnote = 0;
if latlim(1) > latlim(2)
   wid = sprintf('%s:%s:latlimReversed',getcomp,mfilename);
   warning(wid, 'First element of latlim must be less than second')
   errnote = 1;
end
if lonlim(1) > lonlim(2)
   wid = sprintf('%s:%s:lonlimReversed',getcomp,mfilename);
   warning(wid, 'First element of lonlim must be less than second')
   errnote = 1;
end
if errnote
   error(['map:' mfilename ':mapformatsError'], 'Check limits')
end

tolerance = 0;
linebreak = sprintf('\n');

if (latlim(1)>maplatlim(2)+tolerance || ...
    latlim(2)<maplatlim(1)-tolerance || ...
    lonlim(1)>maplonlim(2)+tolerance || ...
    lonlim(2)<maplonlim(1)-tolerance)
    wid = sprintf('%s:%s:limitsOutsideDataset',getcomp,mfilename);
    warning(wid, [ ...
        'Requested latitude or longitude limits are off the map' linebreak ...
        ' latlim for this dataset is ' ...
        mat2str( [maplatlim(1) maplatlim(2)],3) linebreak ...
        ' lonlim for this dataset is '...
        mat2str( [maplonlim(1) maplonlim(2)],3) ...
        ])
    map=[];
    refvec = [];
    return
end

% warn = 0;
if latlim(1) < maplatlim(1)-tolerance
    latlim(1) = maplatlim(1);  % warn = 1;
end

if latlim(2) > maplatlim(2)+tolerance
    latlim(2) = maplatlim(2);  % warn = 1;
end

if lonlim(1) < maplonlim(1)-tolerance
    lonlim(1) = maplonlim(1);  % warn = 1;
end

if lonlim(2) > maplonlim(2)+tolerance
    lonlim(2) = maplonlim(2);  % warn = 1;
end

% if warn
%    warning([ ...
%          'Requested latitude or longitude limits exceed map limits' linebreak ...
%          ' latlim for this dataset is ' ...
%          mat2str( [maplatlim(1) maplatlim(2)],3) linebreak ...
%          ' lonlim for this dataset is '...
%          mat2str( [maplonlim(1) maplonlim(2)],3) ...
%       ])
% end


% convert lat and lonlim to column and row indices
% DTED used to do this:
%   [clim,rlim] = yx2rc(latlim,lonlim,lato,lono,dlat0,dlon0);
% which was equivalent to:
%   clim = ceil( 1 + (latlim - lato)/dlat0 );
%   rlim = ceil( 1 + (lonlim - lono)/dlon0 );
% But it's clear that we need to "snap down" with floor for the lower
% limits, rather than "snap up" with ceil:

clim = [floor(1.5 + (latlim(1) - lato)/dlat0) ...
         ceil(0.5 + (latlim(2) - lato)/dlat0)];
     
rlim = [floor(1.5 + (lonlim(1) - lono)/dlon0) ...
         ceil(0.5 + (lonlim(2) - lono)/dlon0)];

readrows = rlim(1):scalefactor:rlim(2);
readcols = clim(1):scalefactor*skipfactor:clim(2);

% Read the elevation grid
map = readgrid(filename,nrows,ncols,readrows,readcols);

% Construct the referencing vector:  Add a half a cell offset to
% account for the difference between elevation profiles, the paradigm in
% the DTED files, and grid cells with values at the middle, which is the
% model for regular data grids. This will lead to the grid extending a
% half a cell outside the nominal data limits.

% In-line version of "rc2yx"
readlat = (readcols-1)*dlat0 + lato;
readlon = (readrows-1)*dlon0 + lono;

refvec = [1/(dlat*scalefactor),...
             max(readlat) + dlat * (scalefactor - 1/2), ...
             min(readlon) - dlon/2];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function z = readgrid(filename, nrows, ncols, readrows, readcols)

% Read the elevation grid, checking to see if there's any padding at the
% end of the file.  If there is, we need to let READMTX know to expect
% it by setting nFileTrailBytes to a positive (and nonzero) value,
% otherwise READMTX will error.  Note that the data are read in
% transposed form by READMTX, so nrows indicates the number of
% "longitude lines" and ncols indicates the number of "latitude points."
precision = 'int16';
machineformat = 'ieee-be';
nheadbytes = 3428;
nRowHeadBytes  = 8;
nRowTrailBytes = 4;
dirlisting = dir(filename);
expectedFileBytes ...
    = nheadbytes + nrows * (nRowHeadBytes + 2*ncols + nRowTrailBytes);
nFileTrailBytes = max(0, dirlisting.bytes - expectedFileBytes);

z = readmtx(filename,nrows,ncols,precision,readrows,readcols, ...
        machineformat,nheadbytes,nRowHeadBytes,nRowTrailBytes, ...
        nFileTrailBytes);

% Transpose the data, and patch it up, if necessary.
z = correctdata(z');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     
function z = correctdata(z)

% Correct for sign-bit encoded data
% ---------------------------------
% According the DTED Specification (MIL-PRF-89020B, Section 3.11.1),
%
%   "Elevations are two-byte integers, high order first, and negatives are
%   signed magnitude. Users may have to swap the bytes and/or convert
%   negatives to the complement they use.  This can be done by putting the
%   low order byte first, then turning of bit 15 (the high order bit), then
%   multiplying by -1.  For positive numbers, only the bytes are switched."
%
% In other words, standard twos-complement encoding is not used for the
% elevation data.  Instead, each negative value is encoded with exactly the
% same bit pattern as the corresponding positive value, except that the
% most significant bit is flipped from a zero to a one.  We can call this
% "sign-bit encoding."
%
% However, in the call to READMTX (and subsequent calls to FREAD), the
% values are decoded into MATLAB arrays as if they were coming from storage
% in twos-complement form.  Thus corrections are needed following the file
% read process.
%
% The following table shows how various elevation data values (left column)
% are stored using sign-bit encoding (middle column).  From left to right
% we show the two most significant bits and then the three least
% significant bits.  The right column shows the value that FREAD returns
% for each bit pattern when told to read signed 16-bit integers (int16).
%
%         DATA VALUE    ENCODED BITS    VALUE FROM FREAD
%         ----------    ------------    ----------------
%              32767        01...111               32767
%              32766        01...110               32766
%              32765        01...101               32765
%                  3        00...011                   3
%                  2        00...010                   2
%                  1        00...001                   2
%                  0        00...000                   0
%                 -1        10...001              -32767
%                 -2        10...010              -32766
%                 -3        10...011              -32765
%                ...             ...                 ...
%             -32765        11...101                  -3
%             -32766        11...110                  -2
%             -32767        11...111                  -1
%
% The differences in negative values indicate the following corection to
% map a negative value obtained from FREAD back to an actual elevation:
%
%   negative_data_value = -((negative_value_from_fread + 32767) + 1)
%
% The arithmetic form here works for computations in int16.  It's good to
% use it so that later we can enhance DTED to return the elevation grid as
% a class int16 array, rather than double, to economize on memory use.

z(z < 0) = -((z(z < 0) + 32767) + 1);

% Check for non-conforming DTED files
% -----------------------------------
% It turns out that we must check for unexpected negative values just in
% case the file writer used twos-complement after all.  There are files
% available via the Internet for which this occurs.  This step is
% conservative:  we do not interpret the data as being in twos complement
% unless the evidence is very strong.  When in doubt, we assume that the
% spec was adhered to.  Here we assume twos complement only after
% determining that all non-null, negative values would otherwise fall below
% -12000 meters.  This should suffice to catch all twos-complement files
% covering any terrain actually present on Earth.
%
% The DTED specification (Section 3.11.2) says that "in practice,
% the terrain elevation values shall not exceed +9,000 meters or -12,000
% meters."  It also indicates (Section 3.11.3.1) that a "null elevation
% value of -32,767 meters is used as a placeholder in the data record."

negatives = (-32767 < z(:) & z(:) < 0);
if any(negatives) && all(z(negatives) < -12000)
    % Assume that this file was written using twos-complement encoding
    % and warn the user about this assumption.
    wid = sprintf('%s:%s:twosComplementDetected',getcomp,mfilename);
    warning(wid,'%s\n%s',...
        'Negative values are being re-mapped assuming twos-complement',...
        'encoding to try to place them in a reasonable range of magnitudes.')
    % Undo the correction applied above.
    z(z < 0) = (-z(z < 0) - 1) - 32767;
end

% Convert the standard DTED null data value to NaN
z(z == -32767) = NaN;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function deg = dmsstr2deg(str)

% Convert a DTED DMS latitude or longitude string to a number in
% ("decimal") degrees.  The string may or may not contain fractions of
% seconds.  It must be possible, for example, to handle both '123456W' and
% '123456.789W'.  In both strings, '12' refers to 12 degrees and '34'
% refers to 34 minutes.  The seconds component is then either 56 or
% 56.789.  Note that it is not possible to use str2angle here because it
% does not allow for fractions of seconds.

t = str2double(str(1:end-1));
D = fix(t/10000);
t = t - 10000 * D;
M = fix(t/100);
S = t - 100 * M;

deg = D + (M + S/60)/60;
if (str(end) == 'S') || (str(end) == 'W'); 
   deg = -deg; 
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function deg = secstr2deg(str)

% Convert a string containing a latitude or longitude interval (offset) in
% tenths of seconds to a number in ("decimal") degrees.

deg = str2double(str)/36000;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function UHLfield = UHLdescription

% User Header Label field(UHL)

UHLfield( 1).length = 3;    UHLfield( 1).name = 'Recognition sentinel';      
UHLfield( 2).length = 1;    UHLfield( 2).name = 'Fixed by standard';     
UHLfield( 3).length = 8;    UHLfield( 3).name = 'Longitude of origin ';      
UHLfield( 4).length = 8;    UHLfield( 4).name = 'Latitude of origin ';       
UHLfield( 5).length = 4;    UHLfield( 5).name = 'Longitude data interval ';      
UHLfield( 6).length = 4;    UHLfield( 6).name = 'Latitude data interval ';       
UHLfield( 7).length = 4;    UHLfield( 7).name = 'Absolute Vertical Accuracy in Meters';      
UHLfield( 8).length = 3;    UHLfield( 8).name = 'Security Code';     
UHLfield( 9).length = 12;   UHLfield( 9).name = 'Unique reference number ';      
UHLfield(10).length = 4;    UHLfield(10).name = 'number of longitude lines ';        
UHLfield(11).length = 4;    UHLfield(11).name = 'number of latitude points ';        
UHLfield(12).length = 1;    UHLfield(12).name = 'Multiple accuracy';     
UHLfield(13).length = 24;   UHLfield(13).name = 'future use';        

for i=1:length(UHLfield);
   UHLfield(i).type = 'char';
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function DSIfield = DSIdescription

% Data Set Identification (DSI) record contents

DSIfield( 1).length = 3;     DSIfield( 1).name = 'Recognition Sentinel';      
DSIfield( 2).length = 1;     DSIfield( 2).name = 'Security Classification Code';      
DSIfield( 3).length = 2;     DSIfield( 3).name = 'Security Control and Release Markings';     
DSIfield( 4).length = 27;    DSIfield( 4).name = 'Security Handling Description';     
DSIfield( 5).length = 26;    DSIfield( 5).name = 'reserved1';     
DSIfield( 6).length = 5;     DSIfield( 6).name = 'DMA series';        
DSIfield( 7).length = 15;    DSIfield( 7).name = 'unique Ref Num';        
DSIfield( 8).length = 8;     DSIfield( 8).name = 'reserved2';     
DSIfield( 9).length = 2;     DSIfield( 9).name = 'Data Edition Number';       
DSIfield(10).length = 1;     DSIfield(10).name = 'Match Merge Version';       
DSIfield(11).length = 4;     DSIfield(11).name = 'Maintenance Date';      
DSIfield(12).length = 4;     DSIfield(12).name = 'Match Merge Date';      
DSIfield(13).length = 4;     DSIfield(13).name = 'Maintenance Description Code';      
DSIfield(14).length = 8;     DSIfield(14).name = 'Producer Code';     
DSIfield(15).length = 16;    DSIfield(15).name = 'reserved3';     
DSIfield(16).length = 9;     DSIfield(16).name = 'Product Specification';     
DSIfield(17).length = 2;     DSIfield(17).name = 'Product Specification Amendment Number';        
DSIfield(18).length = 4;     DSIfield(18).name = 'Date of Product Specification';     
DSIfield(19).length = 3;     DSIfield(19).name = 'Vertical Datum ';       
DSIfield(20).length = 5;     DSIfield(20).name = 'Horizontal Datum Code ';        
DSIfield(21).length = 10;    DSIfield(21).name = 'Digitizing Collection System';      
DSIfield(22).length = 4;     DSIfield(22).name = 'Compilation Date';      
DSIfield(23).length = 22;    DSIfield(23).name = 'reserved4';     
DSIfield(24).length = 9;     DSIfield(24).name = 'Latitude of origin';        
DSIfield(25).length = 10;    DSIfield(25).name = 'Longitude of origin ';      
DSIfield(26).length = 7;     DSIfield(26).name = 'Latitude of SW corner ';        
DSIfield(27).length = 8;     DSIfield(27).name = 'Longitude of SW corner ';       
DSIfield(28).length = 7;     DSIfield(28).name = 'Latitude of NW corner ';        
DSIfield(29).length = 8;     DSIfield(29).name = 'Longitude of NW corner ';       
DSIfield(30).length = 7;     DSIfield(30).name = 'Latitude of NE corner ';        
DSIfield(31).length = 8;     DSIfield(31).name = 'Longitude of NE corner ';       
DSIfield(32).length = 7;     DSIfield(32).name = 'Latitude of SE corner ';        
DSIfield(33).length = 8;     DSIfield(33).name = 'Longitude of SE corner ';       
DSIfield(34).length = 9;     DSIfield(34).name = 'Clockwise orientation angle ';      
DSIfield(35).length = 4;     DSIfield(35).name = 'Latitude interval ';        
DSIfield(36).length = 4;     DSIfield(36).name = 'Longitude interval ';       
DSIfield(37).length = 4;     DSIfield(37).name = 'Number of Latitude lines';      
DSIfield(38).length = 4;     DSIfield(38).name = 'Number of Longitude lines';     
DSIfield(39).length = 2;     DSIfield(39).name = 'Partial Cell Indicator ';       
DSIfield(40).length = 101;   DSIfield(40).name = 'reserved5';     
DSIfield(41).length = 100;   DSIfield(41).name = 'Reserved for producing nation use ';        
DSIfield(42).length = 156;   DSIfield(42).name = 'reserved6';     

for i=1:length(DSIfield);
   DSIfield(i).type = 'char';
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ACCfield = ACCdescription

% Acuraccy field(ACC) record

ACCfield( 1).length = 3;        ACCfield( 1).name = 'Recognition Sentinel';      
ACCfield( 2).length = 4;        ACCfield( 2).name = 'Absolute Horizontal Accuracy ';     
ACCfield( 3).length = 4;        ACCfield( 3).name = 'Absolute Vertical Accuracy ';       
ACCfield( 4).length = 4;        ACCfield( 4).name = 'Relative Horizontal Accuracy';      
ACCfield( 5).length = 4;        ACCfield( 5).name = 'Relative Vertical Accuracy';        
ACCfield( 6).length = 4;        ACCfield( 6).name = 'reserved1';     
ACCfield( 7).length = 1;        ACCfield( 7).name = 'reserved2';     
ACCfield( 8).length = 31;       ACCfield( 8).name = 'reserved3';     
ACCfield( 9).length = 2;        ACCfield( 9).name = 'Multiple Accuracy Outline Flag';        
ACCfield(10).length = 4;        ACCfield(10).name = 'Sub Absolute Horizontal Accuracy ';     
ACCfield(11).length = 4;        ACCfield(11).name = 'Sub Absolute Vertical Accuracy';        
ACCfield(12).length = 4;        ACCfield(12).name = 'Sub Relative Horizontal Accuracy';      
ACCfield(13).length = 4;        ACCfield(13).name = 'Sub Relative Vertical Accuracy';        
ACCfield(14).length = 14*(2+9+10);      ACCfield(14).name = 'Sub Region Outlines';       
ACCfield(15).length = 18;       ACCfield(15).name = 'reserved4';     
ACCfield(16).length = 69;       ACCfield(16).name = 'reserved5';     

for i=1:length(ACCfield);
   ACCfield(i).type = 'char';
end

%--------------------------------------------------------------------------

function [latdir,londir] = dteddirs(uniquelons,uniquelats,ext)

hWestEast = 'wee';
londir{length(uniquelons)} = [];
for i = 1:length(uniquelons)                                              
    londir{i} = sprintf('dted%c%c%03d%c', filesep,...
        hWestEast(2+sign(uniquelons(i))), abs(uniquelons(i)), filesep);                                             
end                                                                            

hSouthNorth = 'snn';
latdir{length(uniquelats)} = [];
for i = 1:length(uniquelats)                                              
    latdir{i} = sprintf( '%c%02d.%s',...
        hSouthNorth(2+sign(uniquelats(i))), abs(uniquelats(i)), ext);                
end

%--------------------------------------------------------------------------

function [map, refvec, UHL, DSI, ACC] = ...
    readFirstNonEmpty(dtedfile, samplefactor, latlim, lonlim)

% get the origin of the first file
[p,latStr] = fileparts(dtedfile{1,1});
[p,lonStr] = fileparts([p '.*']);
lat0 = str2double(latStr(2:end));
lon0 = str2double(lonStr(2:end));
if strcmp(latStr(1),'s') == 1 || strcmp(latStr(1),'S') == 1
    lat0 = -lat0;
end
if strcmp(lonStr(1),'w') == 1 || strcmp(lonStr(1),'W') == 1
    lon0 = -lon0;
end

% read first non-empty file
nonEmptyFileIndx = [];
for k = 1:numel(dtedfile)
    if exist(dtedfile{k},'file') == 2
        nonEmptyFileIndx = k;
        break;
    end
end
if isempty(nonEmptyFileIndx)
    error(['map:' mfilename ':mapformatsError'], ...
        'No data files for requested area.')
end
[map, refvec, UHL, DSI, ACC] = ...
    dted(dtedfile{nonEmptyFileIndx(1)},samplefactor);

% latitude and longitude limits
dlat = secstr2deg(DSI.Latitudeinterval);
dlon = secstr2deg(DSI.Longitudeinterval);

% map limits for first file read
maplatlim = [ dmsstr2deg(DSI.LatitudeofSWcorner)...
              dmsstr2deg(DSI.LatitudeofNWcorner) ];

% adjust the refvec
topLat   = min([latlim(2) lat0 + diff(maplatlim)]);
leftLon  = max([lonlim(1) lon0]);
refvec(2) = topLat  + dlat * (samplefactor - 1/2);
refvec(3) = leftLon - dlon/2;
