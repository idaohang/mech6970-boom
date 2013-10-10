function paste(this,xShift,yShift)

% Copyright 1996-2008 The MathWorks, Inc.
% $Revision: 1.1.6.5 $  $Date: 2008/10/26 14:25:54 $

position = get(this.hText,'Position');

set(this.hText, ...
    'Position', position + [xShift yShift 0], ...
    'Visible','on')
