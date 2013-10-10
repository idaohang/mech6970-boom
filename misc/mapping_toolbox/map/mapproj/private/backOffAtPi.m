function lon = backOffAtPi(lon, epsilon)

% Back off at +/- 180 degrees to ensure that points originally at +180
% are distinguished from points originally at -180 in the inverse
% projection.  This is necessary for the correct working of the "undo
% trimming and clipping" step, which is applied when direction is
% 'inverse'.

% Copyright 2006 The MathWorks, Inc.
% $Revision: 1.1.6.1 $  $Date: 2006/06/15 20:13:42 $

indx = find(abs(pi - abs(lon)) <= epsilon);
if ~isempty(indx)
    lon(indx) = lon(indx) - sign(lon(indx))*epsilon;
end
