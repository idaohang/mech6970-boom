function cometm(lat,lon,p)
%COMETM Project 2-D comet plot on map axes
%
%  COMETM(lat,lon) projects a comet plot in two dimensions
%  (z = 0) on the current map axes.  A comet plot is an animated
%  graph in which a circle (the comet head) traces the data points on the
%  screen.  The comet body is a trailing segment that follows the head.
%  The tail is a solid line that traces the entire function.  The
%  lat and lon vectors must be in the same units as specified in
%  the map structure.
%
%  COMETM(lat,lon,p) uses the input p to specify a comet body size of
%  p*length(lat)
%
%  See also  COMET3M, COMET, COMET3.

% Copyright 1996-2007 The MathWorks, Inc.
% $Revision: 1.8.4.4 $  $Date: 2007/11/09 20:27:01 $
% Written by:  E. Byrns, E. Brown

if nargin == 0
    errorOnMissingUI(mfilename)
end

checknargin(1,3,nargin,mfilename);
if nargin == 2
    p = [];
end

%  Argument tests
if any([min(size(lat))    min(size(lon))   ] ~= 1) || ...
   any([ndims(lat) ndims(lon)] > 2)

   eid = sprintf('%s:%s:nonVector', getcomp, mfilename);
   error(eid,'%s','Data inputs must be vectors')

elseif ~isequal(size(lat),size(lon))
   eid = sprintf('%s:%s:invalidDimensions', getcomp, mfilename);
   error(eid,'%s','Inconsistent dimensions on input data')

elseif length(p) > 1
   eid = sprintf('%s:%s:invalidTailLength', getcomp, mfilename);
   error(eid,'%s','Tail Length input must be a scalar')
end

%  Test for a map axes and get the map structure
mstruct = gcm;

%  Project the line data
[x,y,z,savepts] = mfwdtran(mstruct,lat,lon,[],'line');

%  Display the comet plot
nextmap;
if isempty(p); 
    comet(x,y)
else
    comet(x,y,p)
end
