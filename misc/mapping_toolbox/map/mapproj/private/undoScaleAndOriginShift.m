function [x, y] = undoScaleAndOriginShift(mstruct, x, y)
% Remove origin shift and return to natural scale

% Copyright 2006 The MathWorks, Inc.
% $Revision: 1.1.6.1 $  $Date: 2006/06/15 20:13:49 $

x = (x - mstruct.falseeasting )/(mstruct.scalefactor);
y = (y - mstruct.falsenorthing)/(mstruct.scalefactor);
