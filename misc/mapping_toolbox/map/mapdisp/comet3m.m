function comet3m(lat,lon,z,p)
%COMET3M Project 3-D comet plot on map axes
%
%  COMET3M(lat,lon,z) projects a comet plot in three dimensions
%  on the current map axes.  A comet plot is an animated graph
%  in which a circle (the comet head) traces the data points on the
%  screen.  The comet body is a trailing segment that follows the head.
%  The tail is a solid line that traces the entire function.  The
%  lat and lon vectors must be in the same units as specified in
%  the map structure.
%
%  COMET3M(lat,lon,z,p) uses the input p to specify a comet body
%  size of p*length(lat)
%
%  See also  COMETM, COMET, COMET3.

% Copyright 1996-2007 The MathWorks, Inc.
% $Revision: 1.14.4.4 $  $Date: 2007/11/09 20:26:59 $
% Written by:  E. Byrns, E. Brown

if nargin == 0
    errorOnMissingUI(mfilename)
end

checknargin(1,4,nargin,mfilename);
if nargin == 3
    p = [];
end

%  Test for scalar z data
%  Comet3 won't accept all data in single z plane, so use z(1) as a work-around
if length(z) == 1;
    z = z(ones(size(lat)));
    z(1) = z(1)-1E-6;
end

%  Argument tests
if any([min(size(lat))    min(size(lon))     min(size(z))] ~= 1) ||...
        any([ndims(lat) ndims(lon)  ndims(z)] > 2)
    eid = sprintf('%s:%s:nonVectorInput', getcomp, mfilename);
    error(eid,'%s','Data inputs must be vectors');

elseif ~isequal(size(lat),size(lon),size(z))
    eid = sprintf('%s:%s:invalidDimensions', getcomp, mfilename);
    error(eid,'%s','Inconsistent dimensions on input data');

elseif length(p) > 1
    eid = sprintf('%s:%s:invalidTailLength', getcomp, mfilename);
    error(eid,'%s','Tail Length input must be a scalar');
end

%  Test for a map axes and get the map structure
mstruct = gcm;

%  Project the line data
[x,y,z,savepts] = mfwdtran(mstruct,lat,lon,z,'line');

%  Display the comet plot
nextmap;
if isempty(p);
    comet3(x,y,z)
else
    comet3(x,y,z,p)
end
