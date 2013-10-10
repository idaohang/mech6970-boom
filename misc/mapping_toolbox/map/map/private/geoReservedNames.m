function reservedNames = geoReservedNames
%GEOFIELDNAMES Return the reserved field names of a geographic data structure 
%
%   NAMES = GEOFIELDNAMES Returns the reserved geographic data structure
%   field names in the cell array NAMES.
%
%   See also GEOATTRIBSTRUCT.

%   Copyright 1996-2006 The MathWorks, Inc.
%   $Revision: 1.1.6.2 $  $Date: 2006/06/15 20:10:59 $

% Return the field names that are reserved
reservedNames = {'Geometry', 'X', 'Y', 'Lat', 'Lon', ...
                 'BoundingBox', 'Height', 'INDEX'};

