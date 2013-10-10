function [x, y] = applyScaleAndOriginShift(mstruct, x, y)
% Apply scale factor, false easting, and false northing

% Copyright 2006 The MathWorks, Inc.
% $Revision: 1.1.6.1 $  $Date: 2006/06/15 20:13:41 $

x = x * mstruct.scalefactor + mstruct.falseeasting;
y = y * mstruct.scalefactor + mstruct.falsenorthing;
