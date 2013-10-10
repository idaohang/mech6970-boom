function hndl = plot3m(varargin)
%PLOT3M Project 3-D lines and points on map axes
%
%  PLOT3M(LAT,LON,Z) projects 3-D line objects onto the current map
%  axes. The input latitude and longitude data must be in the same units
%  as specified in the current map axes.  The units of Z are arbitrary,
%  except when using the 'globe' projection.  In the case of 'globe', Z
%  should have the same units as the radius of the earth or semimajor
%  axis specified in the 'geoid' (reference ellipsoid) property of the
%  map axes.  This implies that for a reference ellipsoid vector of
%  [1 0] (a unit sphere), the units of Z are earth radii.  PLOT3M will
%  clear the current map if the hold state is off.
%
%  PLOT3M(LAT,LON,Z,'LineSpec') uses any valid LineSpec string to
%  display the line object.
%
%  PLOT3M(LAT,LON,Z,'PropertyName',PropertyValue,...) uses the line
%  object properties specified to display the line objects.  Except for
%  xdata, ydata and zdata, all line properties, and styles available
%  through PLOT3 are supported by PLOT3M.
%
%  h = PLOT3M(...) returns the handles to the line objects displayed.
%
%  See also PLOTM, PLOT3, LINEM, LINE.

% Copyright 1996-2009 The MathWorks, Inc.
% $Revision: 1.9.4.6 $  $Date: 2009/05/14 17:06:25 $
% Written by:  E. Byrns, E. Brown

if nargin == 0
    linem;
    return
end

error(nargchk(3, inf, nargin, 'struct'))

%  Display the map
nextmap(varargin)
hndl0 = linem(varargin{:});

%  Set handle return argument if necessary
if nargout == 1
    hndl = hndl0;
end
