function [x,y] = grn2eqa(lat,lon,origin,geoid,units)
%GRN2EQA  Convert from Greenwich to equal area coordinates
%
%  [x,y] = GRN2EQA(lat,lon) transforms from Greenwich latitude and
%  longitude coordinates to equal area cylindrical coordinates.
%  The transformation assumes that the origin of the cylindrical
%  coordinate plane is at (lat,lon) = (0,0).  The inputs are in degrees.
%
%  [x,y] = GRN2EQA(lat,lon,origin) assumes that the origin of the
%  cylindrical coordinate plate is given by the input origin.
%  This input is of the form [lat0,lon0] where the 2 elements are
%  in degrees.
%
%  [x,y] = GRN2EQA(lat,lon,origin,geoid) assumes that the data
%  is distributed on the ellipsoid defined by the input geoid.
%  The geoid vector is of the form [semimajor axes, eccentricity].
%  If omitted, the unit sphere, geoid = [1 0], is assumed.
%
%  [x,y] = GRN2EQA(lat,lon,origin,'units') and
%  [x,y] = GRN2EQA(lat,lon,origin,geoid,'units') assumes that the
%  input data lat, lon and origin are in the angle units specified
%  by 'units'.  If omitted, 'degrees' are assumed.
%
%  mat = GRN2EQA(...) returns a single output matrix,
%  where mat = [x y].
%
%  See also EQA2GRN.

% Copyright 1996-2007 The MathWorks, Inc.
% $Revision: 1.9.4.3 $  $Date: 2007/11/09 20:24:04 $
% Written by:  E. Byrns, E. Brown

error(nargchk(2, 5, nargin, 'struct'))

if nargin == 2
    origin = [];   units = [];    geoid = [];
elseif nargin == 3
    units = [];    geoid = [];
elseif nargin == 4 && ischar(geoid)
    units = geoid;  geoid = [];
elseif nargin == 4 && ~ischar(geoid)
    units = [];
end


%  Empty argument tests
if isempty(units)
    units  = 'degrees';
end
if isempty(geoid)
    geoid  = [1 0];
end
if isempty(origin)
    origin = [0 0 0];
end

%  Dimension tests
if ~isequal(size(lat),size(lon))
    error(['map:' mfilename ':mapError'], ...
        'Inconsistent dimensions on latitude and longitude')
end

%  Test the geoid input
geoid = geoidtst(geoid);

%  Transform the location to an equal area coordinate frame
[x0,y0] = eqacalc(lat,lon,origin,'forward',units,geoid);

%  Set the output arguments
if nargout <= 1
    x = [x0 y0];
else
    x = x0;
    y = y0;
end
