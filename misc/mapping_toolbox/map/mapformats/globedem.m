function [Z, refvec] = globedem(varargin)
%GLOBEDEM Read Global Land One-km Base Elevation (GLOBE) data
%
%  [Z, REFVEC] = GLOBEDEM(FILENAME, SAMPLEFACTOR) reads the GLOBE DEM files
%  and returns the result as a regular matrix map.  The FILENAME is given
%  as a string which does not include an extension.  GLOBEDEM first reads
%  the ESRI header file found in the subdirectory '/esri/hdr/' and then the
%  binary data file filename.  If the files are not found on the MATLAB
%  path, they may be selected interactively. SAMPLEFACTOR is an integer,
%  which when equal to 1 gives the data at its full resolution. When
%  SAMPLEFACTOR is an integer n larger than one, every nth point is 
%  returned. The data grid, Z, is returned as an  array of elevations and 
%  associated referencing vector. Elevations are given in meters above mean
%  sea level using WGS 84 as a horizontal datum.
%
%  [Z, REFVEC] = GLOBEDEM(FILENAME, SAMPLEFACTOR, LATLIM, LONLIM) allows a
%  subset of the map data to be read. The limits of the desired data are
%  specified as vectors of latitude and longitude in degrees. The elements
%  of LATLIM and LONLIM  must be in ascending order.
%
%  [Z, REFVEC] = GLOBEDEM(DIRNAME, SAMPLEFACTOR, LATLIM, LONLIM) reads and 
%  concatenates data from multiple files within a GLOBE directory tree. The
%  DIRNAME input is a string with the name of the directory which contains
%  both the uncompressed files data files and the ESRI header files.  
%
%  For details on locating GLOBEDEM data for download over the Internet,
%  see the following documentation at the MathWorks web site:
%
%  <a href="matlab:
%  web('http://www.mathworks.com/support/tech-notes/2100/2101.html#globe')
%  ">http://www.mathworks.com/support/tech-notes/2100/2101.html</a>
%
%  GLOBE DEM files are binary. No line ending conversion should be performed
%  during transfer or decompression. The ESRI header files are ascii text.
%  Line ending conversion can be applied.
%
%  See also GLOBEDEMS, GTOPO30, DTED, SATBATH, TBASE, USGSDEM.

% Copyright 1996-2007 The MathWorks, Inc.
% $Revision: 1.1.6.5 $ $Date: 2007/10/10 20:49:03 $

name = varargin{1};
if ~isempty(name) & exist(name,'dir') == 7
   if nargin < 4
      error(['map:' mfilename ':mapformatsError'], 'Latlim and lonlim required for directory calling form')
   end
  [Z, refvec] = globedemc(varargin{:});
else 
  [Z, refvec] = globedemf(varargin{:});
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [map,maplegend] = globedemf(fname,scalefactor,latlim,lonlim)


% works for 1.0 data file
%

% use with esri hdr files

if nargin < 1; fname = ''; end
if nargin < 2; scalefactor = 20; end
if nargin < 3; latlim = [-90 90]; end
if nargin < 4; lonlim = [-180 180]; end

latlim = latlim(:)';
lonlim = lonlim(:)';

lonlim = npi2pi(lonlim);

% check input arguments

if length(scalefactor) > 1; 		error(['map:' mfilename ':mapformatsError'], 'SAMPLEFACTOR must be a scalar'); end
if ~isequal([1 2],size(latlim));	error(['map:' mfilename ':mapformatsError'], 'latlim must be a two element vector in units of degrees'); end 
if ~isequal([1 2],size(lonlim)); 	error(['map:' mfilename ':mapformatsError'], 'lonlim must be a two element vector in units of degrees'); end 


%  Open ascii header file and read information

filename = [fname '.hdr'];
fid = fopen(filename,'r');
if fid==-1
   
   % try drilling down to esri/hdr subdirectory. This works 
   % if a full filename has been provided
   [thispth,thisfname,thisext] = fileparts(fname);
   filename = fullfile(thispth,'esri','hdr',[thisfname '.hdr']);
   
   fid = fopen(filename,'r');
   
   if fid==-1
      [filename, path] = uigetfile('*.hdr', ['Select the Globe ESRI header file (' thisfname '.hdr)']);
      if filename == 0 ; return; end
      filename = [path filename];
      fid = fopen(filename,'r');
   end
end

nrows = [];
ncols = [];
nodata = NaN;
ulxmap = [];
ulymap = [];
xdim = [];
ydim = [];

eof = 0;
while ~eof
	str = fscanf(fid,'%s',1);
	switch lower(str)
		case 'nrows', nrows = fscanf(fid,'%d',1);
		case 'ncols', ncols = fscanf(fid,'%d',1);
		case 'nodata', nodata = fscanf(fid,'%d',1);
		case 'ulxmap', ulxmap = fscanf(fid,'%f',1);
		case 'ulymap', ulymap = fscanf(fid,'%f',1);
		case 'xdim', xdim = fscanf(fid,'%f',1);
		case 'ydim', ydim = fscanf(fid,'%f',1);
		case '', eof = 1;
		otherwise, fscanf(fid,'%s',1);
	end
end
fclose(fid);

% other information about the file

precision = 'int16';
machineformat = 'ieee-le';

lato = ulymap;
lono = ulxmap;

dlat = -ydim;
dlon = xdim;

% convert lat and lonlim to column and row indices

[clim,rlim] = yx2rc(lonlim(:),latlim(:),lono,lato,dlon,dlat);

% ensure matrix coordinates are within limits

rlim = [max([1,min(rlim)]) min([max(rlim),nrows])];
clim = [max([1,min(clim)]) min([max(clim),ncols])];

rlim = sort(flipud(rlim(:))');

readrows = rlim(1):scalefactor:rlim(2);
readcols = clim(1):scalefactor:clim(2);

readcols = mod(readcols,ncols); readcols(readcols == 0) = ncols;

% extract the map matrix
map = readmtx(fname,nrows,ncols,precision,readrows,readcols,machineformat);
map = flipud(map);
if ~isempty(map); map(map==nodata) = NaN; end

% Construct the map legend. 
[la1,lo1] = rc2yx(rlim,clim,lato,lono,dlat,dlon);

maplegend = [abs(1/(dlat*scalefactor)) la1(1)-dlat/2 lo1(1)-dlon/2 ];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function str = addequalchar(str)

% ADDEQUALCHAR replace the first space in a string with the equal sign

indx = findstr(str,char(32)); % indices of space characters
str(indx(1)) = '=';
str(end+1) = ';';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [map,maplegend] = globedemc(dname,scalefactor,latlim,lonlim)
%GLOBEDEMD read and concatenate GLOBE DEM files from a directory

[fnames,YMIN,YMAX,XMIN,XMAX,rtile,ctile] = textread('globedems.dat','%s %d %d %d %d %d %d');

% case where dateline is not crossed
if lonlim(1) <= lonlim(2)
	do = ...
	 find( ...
			(...
			(latlim(1) <= YMIN & latlim(2) >= YMAX) | ... % tile is completely within region
			(latlim(1) >= YMIN & latlim(2) <= YMAX) | ... % region is completely within tile
			(latlim(1) >  YMIN & latlim(1) <  YMAX) | ... % min of region is on tile
			(latlim(2) >  YMIN & latlim(2) <  YMAX)   ... % max of region is on tile
			) ...
				&...
			(...
			(lonlim(1) <= XMIN & lonlim(2) >= XMAX) | ... % tile is completely within region
			(lonlim(1) >= XMIN & lonlim(2) <= XMAX) | ... % region is completely within tile
			(lonlim(1) >  XMIN & lonlim(1) <  XMAX) | ... % min of region is on tile
			(lonlim(2) >  XMIN & lonlim(2) <  XMAX)   ... % max of region is on tile
			)...
		);
end


% append root directory and check to see if required files exist
for i = 1:length(do)
   ffname{i} = fullfile(dname,fnames{do(i)});
end	

% assume files exist
fileexist = 1;
for i=1:length(do)
   if ~exist(ffname{i},'file')
      warning('map:globedem:fileNotFound','%s does not exist', ffname{i})
      fileexist = 0;
   end 
end

% exit if not all files exist
if ~fileexist
   error(['map:' mfilename ':mapformatsError'], 'GLOBE DEM file not found.')
end

dortiles = unique(rtile(do));
doctiles = unique(ctile(do));

% read and concatenate separate files
k=0;
map = [];

for i=1:length(dortiles)
   
   rowmap = [];
   
   for j = 1:length(doctiles)
      k = k+1;
      fname = fullfile(dname,fnames{do(k)});
      [tilemap,tilemaplegend] = globedem(fname,scalefactor,latlim,lonlim);
      
      rowmap = [rowmap tilemap];
      if k==1
         maplegend = tilemaplegend;
      end
      
   end
   
   map = [rowmap;map];
   
end

   
      
