function deleteDownload(filename)
%DELETEDOWNLOAD Delete a temporary filename from the system.

%   Copyright 1996-2003 The MathWorks, Inc.
%   $Revision: 1.1.10.3 $ $Date: 2005/11/15 01:07:14 $

try
    delete(filename);
catch
    wid = sprintf('%s:%s:warnDelete', getcomp, mfilename);
    s=sprintf('Unable to delete temporary file "%s".', filename);
    warning(wid, s);
end
