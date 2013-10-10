function lat = backOffAtPoles(lat, epsilon)

% Back off of the +/- 90 degree points.  This allows the differentiation
% of longitudes at the poles of the transformed coordinate system.

% Copyright 2006 The MathWorks, Inc.
% $Revision: 1.1.6.1 $  $Date: 2006/06/15 20:13:43 $

indx = find(abs(pi/2 - abs(lat)) <= epsilon);
if ~isempty(indx)
    lat(indx) = (pi/2 - epsilon) * sign(lat(indx));
end
