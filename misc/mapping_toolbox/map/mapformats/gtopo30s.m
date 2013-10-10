function fname= gtopo30s(latlim,lonlim)
%GTOPO30S GTOPO30 data filenames for latitude-longitude quadrangle
%
%  FNAME = GTOPO30S(LATLIM,LONLIM) returns a cellarray of the file names
%  covering the geographic region for GTOPO30 digital elevation maps (also
%  referred to as "30-arc second" DEMs).  The region is specified by scalar
%  latitude and longitude points, or two element vectors of latitude and
%  longitude limits in units of degrees.
%
%  See also GTOPO30.

% Written by:  A. Kim, W. Stumpf, L. Job
% Copyright 1996-2007 The MathWorks, Inc.
% $Revision: 1.1.6.4 $ $Date: 2007/10/10 20:49:07 $

% ensure row vectors
latlim = latlim(:)';
lonlim = lonlim(:)';

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

fid = fopen('gtopo30s.dat','r');
if fid==-1
	error(['map:' mfilename ':mapformatsError'], 'Couldn''t open gtopo30s.dat')
end

% preallocate bounding rectangle data for speed

YMIN = zeros(1,33); YMAX = YMIN;
XMIN = YMIN; XMAX = YMIN;

% read names and bounding rectangle limits

for n=1:33
	fnames{n,1} = fscanf(fid,'%s',1);
	YMIN(n) = fscanf(fid,'%d',1);
	YMAX(n) = fscanf(fid,'%d',1);
	XMIN(n) = fscanf(fid,'%d',1);
	XMAX(n) = fscanf(fid,'%d',1);
end
fclose(fid);

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
	
% case where the dateline is crossed
if lonlim(1) > lonlim(2)
	lmin = lonlim(1); lmax = lonlim(2);
	lonlim(2) = 180;
	% do eastern side of the dateline first
	doEAST = ...
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
	% do western side of the dateline second
	lonlim(1) = -180; lonlim(2) = lmax;	
	doWEST = ...
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
	% concatenate indices
	do = [doEAST doWEST];
end
	
if ~isempty(do)
	fname = fnames(do);
else
	fname = [];
end

