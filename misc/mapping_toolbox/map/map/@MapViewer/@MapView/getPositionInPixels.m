function p = getPositionInPixels(this)

% Copyright 1996-2008 The MathWorks, Inc.
% $Revision: 1.1.6.3 $  $Date: 2008/10/02 18:56:41 $

oldUnits = get(this.Figure,'Units');
set(this.Figure,'Units','pixels')
p = get(this.Figure,'Position');
set(this.Figure,'Units',oldUnits)
