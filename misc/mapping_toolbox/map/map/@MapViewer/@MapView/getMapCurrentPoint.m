function p = getMapCurrentPoint(this)

% Copyright 1996-2008 The MathWorks, Inc.
% $Revision: 1.1.6.3 $  $Date: 2008/10/26 14:26:33 $

p = get(this.getAxes(),'CurrentPoint');
p = [p(1) p(3)];
