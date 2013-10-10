function errorOnMissingUI(callerName)

% Copyright 2006 The MathWorks, Inc.
% $Revision: 1.1.6.1 $  $Date: 2006/10/14 12:23:44 $

eid = sprintf('%s:%s:missingUI', getcomp, callerName);
error(eid,...
    'Calling %s without arguments to open a dialog box is no longer supported.\nIf you want to enter parameters for %s interactively, create your own dialog box.',...
    upper(callerName),upper(callerName))
