function [outlat,outlon]=track(str,lat,lon,in3,in4,in5)
%TRACK  Track segments to connect navigational waypoints
%
%  [lat,lon] = TRACK(lat0,lon0) connects waypoints given in navigation
%  track format with track segments.  The track segments are rhumb
%  lines, which is the navigationally common method.
%
%  [lat,lon] = TRACK(lat0,lon0,geoid) computes the rhumb line tracks
%  on the ellipsoid defined by the input geoid. The geoid vector
%  is of the form [semimajor axes, eccentricity].  If omitted,
%  the unit sphere, geoid = [1 0], is assumed.
%
%  [lat,lon] = TRACK(lat0,lon0,'units') and
%  [lat,lon] = TRACK(lat0,lon0,geoid,'units') are all valid
%  calling forms, which use the inputs 'units' to define the angle
%  units of the inputs and outputs.  If omitted, 'degrees' are assumed.
%
%  [lat,lon] = TRACK(lat0,lon0,geoid,'units',npts) uses the
%  input npts to determine the number of points per track computed.
%  The input npts is a scalar, and if omitted, npts = 30.
%
%  [lat,lon] = TRACK(mat) and [lat,lon] = TRACK(mat,'units') are
%  valid calling forms, where mat = [lat0 lon0].
%
%  [lat,lon] = TRACK('track',...) uses the 'track' string to define
%  either a great circle or rhumb line tracks.  If 'track' = 'gc',
%  then the great circle tracks are computed.  If 'track' = 'rh', then
%  the rhumb line tracks are computed.  If omitted, 'rh' is assumed.
%
%  mat = TRACK(...) returns a single output argument where mat = [lat lon].
%
%  See also GCWAYPTS, TRACK1, TRACK2, TRACKG.

% Copyright 1996-2007 The MathWorks, Inc.
% $Revision: 1.11.4.5 $  $Date: 2007/11/09 20:25:40 $
% Written by:  E. Brown, E. Byrns

if nargin == 0 || (nargin == 1 && ischar(str))
    
	 error(['map:' mfilename ':mapError'], ...
         'Incorrect number of arguments')

elseif (nargin == 1 && ~ischar(str)) || (nargin == 2 && ischar(str))

    if ~ischar(str)       %  Shift inputs since str omitted by user
        lat = str;
        str = [];
    end

    if size(lat,2) ~= 2 || ndims(lat) > 2
        error(['map:' mfilename ':mapError'], ...
            'Input matrix must have two columns [lat lon]')
    else
        lon = lat(:,2);
        lat = lat(:,1);
    end

    geoid = [];
    units = [];
    npts  = [];

elseif (nargin == 2 && ~ischar(str)) || (nargin == 3 && ischar(str))

    if ~ischar(str)       %  Shift inputs since str omitted by user
        lon = lat;
        lat = str;
        str = [];
    end

    if ischar(lon)         %  track(str,mat,'units')  usage
        if size(lat,2) ~= 2 || ndims(lat) > 2
            error(['map:' mfilename ':mapError'], ...
                'Input matrix must have two columns [lat lon]')
        else
            units = lon;
            lon = lat(:,2);
            lat = lat(:,1);
            geoid = [];
            npts  = [];
        end

    else             %  track(str,lat,lon)  usage
        geoid = [];
        units = [];
        npts  = [];
    end

elseif (nargin == 3 && ~ischar(str)) || (nargin == 4 && ischar(str))

    if ~ischar(str)       %  Shift inputs since str omitted by user
        in3  = lon;
        lon = lat;
        lat  = str;
        str  = [];
    end

    if ischar(in3)
        geoid = [];
        units = in3;
        npts  = [];
    else
        geoid = in3;
        units = [];
        npts  = [];
    end

elseif (nargin == 4 && ~ischar(str)) || (nargin == 5 && ischar(str))

    if ~ischar(str)       %  Shift inputs since str omitted by user
        in4  = in3;
        in3  = lon;
        lon  = lat;
        lat  = str;
        str  = [];
    end

    geoid = in3;
    units = in4;
    npts  = [];


elseif (nargin == 5 && ~ischar(str)) || (nargin == 6 && ischar(str))

    if ~ischar(str)       %  Shift inputs since str omitted by user
        in5  = in4;
        in4  = in3;
        in3  = lon;
        lon  = lat;
        lat  = str;
        str  = [];
    end

    geoid = in3;
    units = in4;
    npts  = in5;

end


%  Test the track string

if isempty(str)
    str = 'rh';       %  Default is rhumb line tracks
else
    validstr = {'gc','rh'};
    indx = find(strcmpi(str,validstr));
    if numel(indx) ~= 1
        error(['map:' mfilename ':mapError'], 'Unrecognized track string')
    else
        str = validstr{indx};
    end
end

%  Empty argument tests.  Set defaults

if isempty(geoid)
    geoid = [0 0];
end

if isempty(units)
    units = 'degrees';
end

if isempty(npts)
    npts  = 30;
end

%  Input dimension tests

if ~isequal(size(lat),size(lon))
    error(['map:' mfilename ':mapError'], ...
        'Inconsistent dimensions for lat and lon inputs.')
elseif any([min(size(lat)) min(size(lon))] ~= 1) || ...
       any([ndims(lat) ndims(lon)] > 2)
    error(['map:' mfilename ':mapError'], ...
        'Latitude and longitude inputs must vectors')
elseif max(size(lat)) == 1
    error(['map:' mfilename ':mapError'], ...
        'At least 2 lats and 2 longs are required')
elseif max(size(npts)) ~= 1
    error(['map:' mfilename ':mapError'], ...
        'Scalar npts required')
end

%  Test the geoid parameter

geoid = geoidtst(geoid);

[lat, lon] = toRadians(units, lat, lon);
[outlat, outlon] = doTrack(str, lat, lon, geoid, npts);
[outlat, outlon] = fromRadians(units, outlat, outlon);

%  Set the output argument if necessary
if nargout < 2
    outlat = [outlat outlon];
end
