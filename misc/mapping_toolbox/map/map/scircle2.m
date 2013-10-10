function [latout,lonout] = scircle2(str,lat1,lon1,lat2,lon2,in5,in6,in7)
%SCIRCLE2  Small circles from center and perimeter
% 
%   [lat,lon] = SCIRCLE2(lat1,lon1,lat2,lon2) computes small circles
%   (on a sphere) with centers at the point lat1, lon1 and points on
%   the circles at lat2, lon2.  The inputs can be scalar or column vectors.
% 
%   [lat,lon] = SCIRCLE2(lat1,lon1,lat2,lon2,geoid) computes the small
%   circle on the ellipsoid defined by the input geoid, rather than
%   assuming a sphere. The geoid vector is of the form [semimajor axis,
%   eccentricity].  If geoid = [], a sphere is assumed.
% 
%   [lat,lon] = SCIRCLE2(lat1,lon1,lat2,lon2,'units') and
%   [lat,lon] = SCIRCLE2(lat1,lon1,lat2,lon2,geoid,'units') are all valid
%   calling forms, which use the inputs 'units' to define the angle
%   units of the inputs and outputs.  If omitted, 'degrees' are assumed.
% 
%   [lat,lon] = SCIRCLE2(lat1,lon1,lat2,lon2,geoid,'units',npts) uses
%   the scalar input npts to determine the number of points per track
%   computed.  The default value of npts is 100.
% 
%   [lat,lon] = SCIRCLE2('track',...) uses the 'track' string to define
%   either a great circle or rhumb line radius.  If 'track' = 'gc',
%   then small circles are computed.  If 'track' = 'rh', then
%   the circles with radii of constant rhumb line distance are computed.
%   If omitted, 'gc' is assumed.
% 
%   mat = SCIRCLE2(...) returns a single output argument where
%   mat = [lat lon].  This is useful if only a single circle is computed.
% 
%   Multiple circles can be defined from a single center point by
%   providing scalar lat1, lon1 inputs and column vectors for the points
%   on the circumference, lat2, lon2.
% 
%   See also SCIRCLE1, SCIRCLEG, TRACK2.

% Copyright 1996-2007 The MathWorks, Inc.
% $Revision: 1.9.4.4 $  $Date: 2007/11/09 20:25:24 $
% Written by:  E. Byrns, E. Brown

if nargin == 0
     error(['map:' mfilename ':mapError'], ...
         'Incorrect number of arguments')
elseif (nargin < 4  && ~ischar(str)) || (nargin == 4 && ischar(str))
	 error(['map:' mfilename ':mapError'], ...
         'Incorrect number of arguments')
elseif (nargin == 4 && ~ischar(str)) || (nargin == 5 && ischar(str))

    if ~ischar(str)       %  Shift inputs since str omitted by user
		lon2 = lat2;
        lat2 = lon1;
		lon1 = lat1;
        lat1 = str;
		str  = [];
    end

	geoid = [];    units = [];    npts  = [];

elseif (nargin == 5 && ~ischar(str)) || (nargin == 6 && ischar(str))

    if ~ischar(str)       %  Shift inputs since str omitted by user
	    in5  = lon2;     lon2 = lat2;
		lat2 = lon1;     lon1 = lat1;
		lat1 = str;      str  = [];
    end

    if ischar(in5)
	    geoid = [];     units = in5;     npts  = [];
    else
	    geoid = in5;    units = [];      npts  = [];
    end


elseif (nargin == 6 && ~ischar(str)) || (nargin == 7 && ischar(str))

    if ~ischar(str)       %  Shift inputs since str omitted by user
	    in6  = in5;      in5  = lon2;
		lon2 = lat2;     lat2 = lon1;
		lon1 = lat1;     lat1 = str;
		str  = [];
    end

    geoid = in5;      units = in6;      npts  = [];


elseif (nargin == 7 && ~ischar(str)) || (nargin == 8 && ischar(str))

    if ~ischar(str)       %  Shift inputs since str omitted by user
	    in7  = in6;      in6  = in5;
	    in5  = lon2;     lon2 = lat2;
		lat2 = lon1;     lon1 = lat1;
		lat1 = str;      str  = [];
    end

    geoid = in5;      units = in6;       npts  = in7;

end


%  Test the track string

if isempty(str)
    str = 'gc';
else
    validstr = {'gc','rh'};
	indx     = strmatch(lower(str),validstr);
    if length(indx) ~= 1
        error(['map:' mfilename ':mapError'], 'Unrecognized track string')
    else
        str = validstr{indx};
    end
end

%  Allow for scalar starting point, but vectorized azimuths.  Multiple
%  circles starting from the same point

if length(lat1) == 1 && length(lon1) == 1 && ~isempty(lat2)
    lat1 = lat1(ones(size(lat2)));   lon1 = lon1(ones(size(lat2)));
end

%  Empty argument tests.  Set defaults

if isempty(geoid);   geoid = [0 0];       end
if isempty(units);   units = 'degrees';   end
if isempty(npts);    npts  = 100;         end

%  Dimension tests

if ~isequal(size(lat1),size(lon1),size(lat2),size(lon2))
      error(['map:' mfilename ':mapError'], ...
          'Inconsistent dimensions on latitude and longitude inputs')
elseif max(size(npts)) ~= 1
       error(['map:' mfilename ':mapError'], 'Scalar npts required')
end

%  Test the geoid parameter

geoid = geoidtst(geoid);

% Ensure that inputs are column vectors

lat1 = lat1(:);    lon1 = lon1(:);
lat2 = lat2(:);    lon2 = lon2(:);

%  Angle unit conversion

[lat1, lon1, lat2, lon2] ...
    = toRadians(units, lat1, lon1, lat2, lon2);

%  Compute azimuth and range

rng = distance(str,lat1,lon1,lat2,lon2,geoid,'radians');
az  = 2*pi;             az = az(ones([size(lat1,1) 1]));

%  Compute circles

[latc,lonc] = scircle1(str,lat1,lon1,rng,az,geoid,'radians',npts);

%  Convert the results to the desired units

[latc, lonc] = fromRadians(units, latc, lonc);

%  Set the output arguments

if nargout <= 1
     latout = [latc lonc];
elseif nargout == 2
     latout = latc;
     lonout = lonc;
end
