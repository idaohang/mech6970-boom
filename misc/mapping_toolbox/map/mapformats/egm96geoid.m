function [map,maplegend] = egm96geoid(scalefactor,latlim,lonlim)
% EGM96GEOID Read 15-minute gridded geoid heights from EGM96
% 
%  [N, REFVEC] = EGM96GEOID(SAMPLEFACTOR) imports global geoid height in
%  meters from the EGM96 geoid model. The data set is gridded at 15-minute
%  intervals, but may be downsampled as specified by the positive integer
%  SAMPLEFACTOR. The result is returned in the regular data grid N along
%  with referencing vector REFVEC. At full resolution (a SAMPLEFACTOR of 1), 
%  N will be 721-by-1441.
%
%  The gridded EGM96 data set must be on your path in a file named
%  'WW15MGH.GRD'. 
% 
%  [N, REFVEC] = EGM96GEOID(SAMPLEFACTOR, LATLIM, LONLIM) imports data for
%  the part of the world within the specified latitude and longitude
%  limits. The limits must be two-element vectors in units of degrees.
%  Longitude limits can be defined in the range [-180 180] or [0 360]. For
%  example, lonlim = [170 190] returns data centered on the dateline, while
%  lonlim = [-10 10] returns data centered on the prime meridian.
%
%  For details on locating EGM96 geoid data for download over the Internet,
%  see the following documentation at the MathWorks web site: 
%
%  <a href="matlab: 
%  web('http://www.mathworks.com/support/tech-notes/2100/2101.html#egm96') 
%  ">http://www.mathworks.com/support/tech-notes/2100/2101.html</a>

% Copyright 1996-2007 The MathWorks, Inc.
% $Revision: 1.1.6.5 $  $Date: 2007/10/10 20:48:59 $
% Written by:  W. Stumpf

if nargin == 1
	latlim = [-90 90]; % degrees
	lonlim = [0 360];
elseif nargin == 3
	latlim = sort(latlim(:)');
	lonlim = npi2pi(lonlim(:)'); % No effort made (yet) to work across the dateline
	if lonlim(2) < lonlim(1); lonlim(2) = lonlim(2) + 360; end
else
    error('map:egm96geoid:invalidArgCount', ...
          'Incorrect number of arguments.')
end

latlim = latlim(:)';
lonlim = lonlim(:)';

% check input arguments

if length(scalefactor) > 1
    error('map:egm96geoid:invalidSamplefactor','SAMPLEFACTOR must be a scalar')
end
if ~isequal([1 2],size(latlim))
    error('map:egm96geoid:invalidLatlim',...
        'latlim must be a two element vector in units of degrees')
end 
if ~isequal([1 2],size(lonlim))
    error('map:egm96geoid:invalidLonlim',...
        'lonlim must be a two element vector in units of degrees')
end 

% Open the file, read the file header, and close it again

filename = 'WW15MGH.GRD';
fid = fopen(filename,'r');

if fid == -1
	filename = lower(filename);
	fid = fopen(filename,'r');
	if fid == -1
		[filename,filepath] = uigetfile(filename,['Where is ',filename,'?']);
		if filename == 0 ; return; end
		fid = fopen([filepath,filename],'r');
		filename = [filepath,filename];
	end
end


fields(1).name = 'south' ; fields(1).type = '%13g'; fields(1).length = 1;
fields(2).name = 'north' ; fields(2).type = '%12g'; fields(2).length = 1;
fields(3).name = 'west'  ; fields(3).type = '%12g'; fields(3).length = 1;
fields(4).name = 'east'  ; fields(4).type = '%12g'; fields(4).length = 1;
fields(5).name = 'dphi'  ; fields(5).type = '%12g'; fields(5).length = 1;
fields(6).name = 'dlam'  ; fields(6).type = '%12g'; fields(6).length = 1;

s = readfields(filename,fields,1,'native',fid);

% find end of file

fseek(fid,0,1);
fsize = ftell(fid);

% close file 

fclose(fid);

% parameters describing file format and contents

lato =   s.north;
lono =   s.west;
dlat = - s.dphi;
dlon =   s.dlam;

nrows = 721;
ncols = 1441;

% file has line breaks within logical records. Compute the number
% of bytes per record including line breaks. We bother with this because
% we can read just the data we want using READMTX.

% for PC line endings, change lineend to 2, and assume 
% last carriage return is followed by a linefeed

if fsize == 9618935 % one line ending character (mac and unix)

	bigblocks = 9;
	lineperbigblock = 20;
	valuesperline = 8;
	blockchars = 73; % bytes
	lineend = 1; % bytes
	lastpoint = 10; % bytes
	nFileTrailBytes = 0;
elseif fsize == 9756647 % two line ending characters (PC)

	bigblocks = 9;
	lineperbigblock = 20;
	valuesperline = 8;
	blockchars = 73; % bytes
	lineend = 2; % bytes
	lastpoint = 10; % bytes
	nFileTrailBytes = 1; % no linefeed
else
	error('map:egm96geoid:invalidFileSize', '%s', ...
        'File size doesn''t match original. Try downloading and \nextracting the ''WW15MGH.GRD'' file as binary.')
end
	
bytesperlat = bigblocks*( lineperbigblock*(blockchars+lineend)  + lineend ) ...
				+ lastpoint + 2*lineend;


% convert lat and lonlim to column and row indices

[clim,rlim] = yx2rc(lonlim(:),latlim(:),lono,lato,dlon,dlat);

% ensure matrix coordinates are within limits

rlim = sort(flipud(rlim(:))');

readrows = rlim(1):scalefactor:rlim(2);
readcols = clim(1):scalefactor:clim(2);

% rewrap areas straddling ends of matrix

readcols = mod(readcols,ncols); readcols(readcols == 0) = ncols;

% take into account repeated column at the ends of the matrix

indx = find(diff(readcols)<1);
if ~isempty(indx);
	readcols(indx(1)+1:end) = readcols(indx(1)+1:end)+1;
	indx2 = find(readcols(indx(1)+1:end)>clim(2));
    readcols(indx(1)+indx2) = [];
end

% extract the map matrix

map = readmtx(filename,nrows,ncols,'%9g',readrows,readcols,'native',74,0,0,nFileTrailBytes,bytesperlat);
map = flipud(map);

% Construct the map legend. Add a half a cell offset to the map legend to 
% account for the difference between lat and long coincident with the value, 
% and matrix cells with values at the middle, which is the model 
% for regular matrix maps. This may lead to the map extending a half a cell
% outside the requested map limits.

[la1,lo1] = rc2yx(rlim,clim,lato,lono,dlat,dlon);

maplegend = [abs(1/(dlat*scalefactor)) la1(1)-dlat/2 lo1(1)-dlon/2 ];

