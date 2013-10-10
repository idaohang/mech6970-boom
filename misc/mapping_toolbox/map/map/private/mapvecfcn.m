function fcn = mapvecfcn(geometry, fcnname)
%MAPVECFCN Get map vector display function for given geometry
%
%   FCN = MAPVECVCN(GEOMETRY, FCNNAME) returns a function handle in FCN
%   based on the GEOMETRY string. Valid GEOMETRY strings are 'point',
%   'multipoint', 'line', or 'polygon'. If the geometry is not determined,
%   an error message is constructed using the name of the calling function,
%   FCNNAME.
%
%   See also MAPSHOW, MAPSTRUCTFCN, MAPSTRUCTSHOW, MAPVECSHOW.

% Copyright 2006 The MathWorks, Inc.
% $Revision: 1.1.6.1 $  $Date: 2006/05/24 03:35:46 $

switch (lower(geometry))

   case {'point', 'multipoint'}
      fcn = @mappoint;

   case 'line'
      fcn = @mapline;

   case 'polygon'
      fcn = @mappolygon;

   otherwise
      eid = sprintf('%s:%s:invalidGeometry', getcomp, fcnname);
      msg = sprintf('%s%s%s\n%s\n%s%s%s', 'Function ', upper(fcnname), ...
         ' expected ''DisplayType'' to be ', ...'
         '''Point'', ''MultiPoint'', ''Line'', or ''Polygon'';', ...
         'instead it was ''', geometry, '''.');
      error(eid,'%s',msg)

end
