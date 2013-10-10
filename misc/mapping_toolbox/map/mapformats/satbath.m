function [latgrat,longrat,mat] = satbath(varargin)
%SATBATH Read 2-minute global terrain/bathymetry from Smith and Sandwell
%
%  [LATGRAT, LONGRAT, Z] = SATBATH reads the global topography file for the
%  entire world (topo_8.2.img), returning every 50th point.  The result is
%  returned as a geolocated data grid. If using a different version of the
%  global topography file, rename it to topo_8.2.img. If the file is not
%  found on the MATLAB path, a dialog will open to request the file.
%
%  [LATGRAT, LONGRAT, Z] = SATBATH(SAMPLEFACTOR) returns the data for the
%  entire world, subsampled by the integer SAMPLEFACTOR.  A samplefactor of
%  n returns every nth point.  The data grid at full resolution has 6336 by
%  10800 points.
%
%  [LATGRAT, LONGRAT, Z] = SATBATH(SAMPLEFACTOR, LATLIM, LONLIM) returns
%  data for the specified region.  The returned data will extend slightly
%  beyond the requested area.  If omitted, the entire area covered by the
%  data file is returned.  The limits are ascending two element vectors in
%  units of degrees.
%
%  [LATGRAT, LONGRAT, Z] = SATBATH(SAMPLEFACTOR, LATLIM, LONLIM, GSIZE)
%  controls the size of the graticule matrices.  GSIZE is a two-element
%  vector containing the number of rows and columns desired.  If omitted, a
%  graticule the size of the map is returned.
%
%  For details on locating SATBATH data for download over the Internet, see
%  the following documentation at the MathWorks web site:
%
%  <a href="matlab:
%  web('http://www.mathworks.com/support/tech-notes/2100/2101.html#scripps')
%  ">http://www.mathworks.com/support/tech-notes/2100/2101.html</a>
%
%  See also TBASE, GTOPO30, EGM96GEOID.

% [nrows,ncols] = satbath('size',SAMPLEFACTOR,LATLIM,LONLIM) returns the
% size of the matrix without extracting the data.

% Copyright 1996-2007 The MathWorks, Inc.
% $Revision $  $Date: 2007/10/10 20:49:12 $

% Calling form to compute size of matrix without extracting data.
request = 'data';
nargi= nargin;
if nargin >= 1 && ischar(varargin{1})
   request = varargin{1};
   varargin(1) = [];
   nargi = nargin-1;
end

switch request
   case{'data','size'}
   otherwise
      error('map:satbath:invalidRequest',...
          'Valid requests are ''data'' and ''size''')
end

if nargi < 1
   scalefactor = 50;
else
   scalefactor = varargin{1};
end

if nargi < 2
   latlim = [-72 72];
else
   latlim = varargin{2};
end 

if nargi < 3
   lonlim = [0 360];
else
   lonlim = varargin{3};
end 

if nargi < 4
   gsize = [];
else
   gsize = varargin{4};
end

% Get lonlim into the correct range.
latlim = latlim(:)';
lonlim = lonlim(:)';

if ~isequal(lonlim,[0 360])
   lonlim = zero22pi(lonlim);
   if lonlim(2) < lonlim(1);
      lonlim(2) = lonlim(2) + 360;
   end
end

% The size of the matrix.
nrows = 6336;
ncols = 10800;

% Map lat and long limits to projected coordinates.
[rlim,clim] = ll2rc(latlim,lonlim);

% Trim latitudes to the limits of the map.
rlim = [max([1,min(rlim)]) min([max(rlim),nrows])];

% Extract the map matrix.
readrows = rlim(1):scalefactor:rlim(2);
readcols = clim(1):scalefactor:clim(2);

% Bail if all we wanted was the size.
switch request
   case 'size'
      [latgrat,longrat] = deal(length(readrows),length(readcols));
      return
end

% Assign filename to the current version of the topography file.
filename = 'topo_8.2.img';

% Verify the file is on the MATLAB path.
if exist(filename,'file') ~= 2
   if exist('topo_6.2.img','file') == 2
      filename = 'topo_6.2.img';
   else
      [filename,filepath] = uigetfile(filename,['Where is ',filename,'?']);
      if filename == 0
         latgrat = [];
         longrat = [];
         mat = [];
         return
      else
         filename = [filepath,filename];
      end
   end
end

% Construct a graticule of row and column indices.
if nargi < 4 		% size(grat) = size(mat)
   [rIndGrat,cIndGrat] = meshgrat(readrows,readcols);
else  	         % texture map the data to a smaller graticule
   [rIndGrat,cIndGrat] = meshgrat([min(readrows) max(readrows)], ...
      [min(readcols) max(readcols)],gsize);
end

% Wrap data that extends beyond the right edge of the map.
readcols = mod(readcols,ncols);
readcols(readcols == 0) = ncols;

% Check for requests straddling the edge of the data.
indx = find(diff(readcols) < 1);

% Read the data.
if isempty(indx)   % no straddle
   mat = readmtx(filename,nrows,ncols,'integer*2',readrows, ...
      readcols,'ieee-be');
else
   mat1 = readmtx(filename,nrows,ncols,'integer*2',readrows, ...
      readcols(1:indx(1)),'ieee-be');
   mat2 = readmtx(filename,nrows,ncols,'integer*2',readrows, ...
      readcols(indx(1)+1:end),'ieee-be');
   mat = [mat1 mat2];
end

% Map row and column graticule to lat and long.
[latgrat,longrat] = rc2ll(rIndGrat,cIndGrat);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [i,j] = ll2rc(lat,lon)
%LL2RC converts lat-long to row-column for the global satellite bathymetry

% a translated excerpt of Davis Sandwell's copyrighted code in
% ftp://topex.ucsd.edu/pub/global_topo_2min/src_img_sun/img2xyt.f

d2r =.0174533;

nrows=6336;
lat1=-72.006;
lon1=0.;
dlon=2/60.;

arg1=log(tan(d2r*(45.+lat1/2)));
arg2=log(tan(d2r*(45.+lat/2)));
i=ceil(nrows+1-(arg2-arg1)/(dlon*d2r));
j=floor((lon-lon1)/dlon+1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [lat,lon] = rc2ll(i,j)
%RC2ll converts lat-long to row-column for the global satellite bathymetry

% a translated excerpt of Davis Sandwell's copyrighted code in
% ftp://topex.ucsd.edu/pub/global_topo_2min/src_img_sun/img2xyt.f

d2r =.0174533;
nrows=6336;
lat1=-72.006;
lon1=0.;
dlon=2/60.;

arg1=d2r*dlon*(nrows-i+.5);
arg2=log(tan(d2r*(45.+lat1/2.)));
term=exp(arg1+arg2);
lat=2.*atan(term)/d2r-90.;
lon=lon1+dlon*(j-.5);
