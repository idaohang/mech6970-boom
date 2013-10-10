function [x, y] = removeDuplicateVertices(x, y)
% Given a pair of vertex arrays, remove any vertex which is nearly
% identical to its successor. That is, if x(k) == x(k+1) and y(k) ==
% y(k+1), or if they match to within 10*eps(x(k)), remove the k-th
% vertex.

% Copyright 2008 The MathWorks, Inc.
% $Revision: 1.1.6.1 $  $Date: 2008/11/24 15:00:19 $

if numel(x) > 1
    epsX = eps(x(1:(end-1)));
    epsY = eps(y(1:(end-1)));
    duplicates = [...    
        (abs(diff(x(:))) < 10*epsX(:)) & ...
        (abs(diff(y(:))) < 10*epsY(:)); false];
    x(duplicates) = [];
    y(duplicates) = [];
end
