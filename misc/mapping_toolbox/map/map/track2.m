function [latout,lonout] = track2(str,lat1,lon1,lat2,lon2,in5,in6,in7)
%TRACK2  Geographic tracks from starting and ending points
% 
%   [lat,lon] = TRACK2(lat1,lon1,lat2,lon2) computes great circle
%   tracks on a sphere starting at the point lat1, lon1 and ending at
%   lat2, lon2.  The inputs can be scalar or column vectors.
% 
%   [lat,lon] = TRACK2(lat1,lon1,lat2,lon2,geoid) computes the great circle
%   track on the ellipsoid defined by the input geoid. The geoid vector
%   is of the form [semimajor axis, eccentricity].  If geoid = [], a sphere
%   is assumed.
% 
%   [lat,lon] = TRACK2(lat1,lon1,lat2,lon2,'units') and
%   [lat,lon] = TRACK2(lat1,lon1,lat2,lon2,geoid,'units') are all valid
%   calling forms, which use the inputs 'units' to define the angle
%   units of the inputs and outputs.  If omitted, 'degrees' are assumed.
% 
%   [lat,lon] = TRACK2(lat1,lon1,lat2,lon2,geoid,'units',npts) uses the
%   scalar input npts to determine the number of points per track computed.
%   The default value of npts is 100.
% 
%   [lat,lon] = TRACK2('track',...) uses the 'track' string to define
%   either a great circle or rhumb line track.  If 'track' = 'gc',
%   then the great circle tracks are computed.  If 'track' = 'rh', then
%   the rhumb line tracks are computed.  If omitted, 'gc' is assumed.
% 
%   mat = TRACK2(...) returns a single output argument where
%   mat = [lat lon].  This is useful if only a single track is computed.
% 
%   Multiple tracks can be defined from a single starting point by
%   providing scalar lat1, lon1 inputs and column vectors for lat2, lon2.
% 
%   See also TRACK1, TRACKG, SCIRCLE2.

% Copyright 1996-2007 The MathWorks, Inc.
% $Revision: 1.9.4.5 $  $Date: 2007/11/09 20:25:44 $
% Written by:  E. Byrns, E. Brown

if nargin == 0
    error(['map:' mfilename ':mapError'], 'Incorrect number of arguments')

elseif (nargin < 4  && ~ischar(str)) || (nargin == 4 && ischar(str))
    error(['map:' mfilename ':mapError'], 'Incorrect number of arguments')

elseif (nargin == 4 && ~ischar(str)) || (nargin == 5 && ischar(str))

    if ~ischar(str)       %  Shift inputs since str omitted by user
        lon2 = lat2;
        lat2 = lon1;
        lon1 = lat1;
        lat1 = str;
        str  = [];
    end

    geoid = [];
    units = [];
    npts  = [];

elseif (nargin == 5 && ~ischar(str)) || (nargin == 6 && ischar(str))

    if ~ischar(str)       %  Shift inputs since str omitted by user
        in5  = lon2;
        lon2 = lat2;
        lat2 = lon1;
        lon1 = lat1;
        lat1 = str;
        str  = [];
    end

    if ischar(in5)
        geoid = [];
        units = in5;
        npts  = [];
    else
        geoid = in5;
        units = [];
        npts  = [];
    end

elseif (nargin == 6 && ~ischar(str)) || (nargin == 7 && ischar(str))

    if ~ischar(str)       %  Shift inputs since str omitted by user
        in6  = in5;
        in5  = lon2;
        lon2 = lat2;
        lat2 = lon1;
        lon1 = lat1;
        lat1 = str;
        str  = [];
    end

    geoid = in5;
    units = in6;
    npts  = [];

elseif (nargin == 7 && ~ischar(str)) || (nargin == 8 && ischar(str))

    if ~ischar(str)       %  Shift inputs since str omitted by user
        in7  = in6;
        in6  = in5;
        in5  = lon2;
        lon2 = lat2;
        lat2 = lon1;
        lon1 = lat1;
        lat1 = str;
        str  = [];
    end

    geoid = in5;
    units = in6;
    npts  = in7;

end


%  Test the track string

if isempty(str)
    str = 'gc';
else
    validstr = {'gc','rh'};
    indx = strmatch(lower(str),validstr);
    if numel(indx) ~= 1
        error(['map:' mfilename ':mapError'], 'Unrecognized track string')
    else
        str = validstr{indx};
    end
end

%  Allow for scalar starting point, but vectorized azimuths.  Multiple
%  tracks starting from the same point

if (numel(lat1) == 1) && (numel(lon1)) == 1 && ~isempty(lat2)
    lat1 = lat1(ones(size(lat2)));
    lon1 = lon1(ones(size(lat2)));
end

%  Empty argument tests.  Set defaults

if isempty(geoid)
    geoid = [0 0];
end

if isempty(units)
    units = 'degrees';
end

if isempty(npts)
    npts  = 100;
end

%  Dimension tests
if ~isequal(size(lat1),size(lon1),size(lat2),size(lon2))
    error(['map:' mfilename ':mapError'], ...
        'Inconsistent dimensions on latitude and longitude inputs')
elseif max(size(npts)) ~= 1
    error(['map:' mfilename ':mapError'], ...
        'Scalar npts required')
end

%  Test the geoid parameter
geoid = geoidtst(geoid);

%  Angle unit conversion
[lat1, lon1, lat2, lon2] = toRadians(units, lat1, lon1, lat2, lon2);

[lattrk, lontrk] = doTrack2(str, lat1, lon1, lat2, lon2, geoid, npts);

%  Convert the results to the desired units
[lattrk, lontrk] = fromRadians(units, lattrk, lontrk);

%  Set the output arguments
if nargout <= 1
     latout = [lattrk lontrk];
elseif nargout == 2
     latout = lattrk;
     lonout = lontrk;
end
