function [fname,qname] = usgsdems(latlim,lonlim)
%USGSDEMS USGS 1-Degree DEM filenames for latitude-longitude quadrangle
%
%   [FNAME, QNAME] = USGSDEMS(LATLIM, LONLIM) returns cellarrays of the file
%   names and quadrangle names covering the geographic region for 1-degree
%   USGS digital elevation maps (also referred to as "3-arc second" or
%   "1:250,000 scale" DEMs).  The region is specified by scalar latitude
%   and longitude points, or two element vectors of latitude and longitude
%   limits in units of degrees.
%
%   See also DEMDATAUI, USGS24KDEM, USGSDEM.

% Copyright 1996-2007 The MathWorks, Inc.
% $Revision: 1.1.6.4 $
% Written by:  A. Kim, W. Stumpf

if nargin~=2
	error(['map:' mfilename ':mapformatsError'], 'Incorrect number of arguments')
end

if  isequal(size(latlim),[1 1])
	latlim = latlim*[1 1];
elseif ~isequal(size(latlim),[1 2])
    error(['map:' mfilename ':mapformatsError'], 'Latitude limit input must be a scalar or 2 element vector')
end

if isequal(sort(size(lonlim)),[1 1])
	lonlim = lonlim*[1 1];
elseif ~isequal(sort(size(lonlim)),[1 2])
    error(['map:' mfilename ':mapformatsError'], 'Longitude limit input must be a scalar or 2 element vector')
end

fid = fopen('usgsdems.dat','r');
if fid==-1
	error(['map:' mfilename ':mapformatsError'], 'Couldn''t open usgsdems.dat')
end

% preallocate bounding rectangle data for speed

YMIN = zeros(1,924); YMAX = YMIN;
XMIN = YMIN; XMAX = YMIN;

% read names and bounding rectangle limits

for n=1:924
	fnames{n,1} = fscanf(fid,'%s',1);
	YMIN(n) = fscanf(fid,'%d',1);
	YMAX(n) = fscanf(fid,'%d',1);
	XMIN(n) = fscanf(fid,'%d',1);
	XMAX(n) = fscanf(fid,'%d',1);
	qnames{n,1} = fscanf(fid,'%s',1);
end
fclose(fid);


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

if ~isempty(do)
	fname = fnames(do);
	qname = qnames(do);
else
	fname = [];
	qname = [];
end
