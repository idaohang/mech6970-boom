function closeFiles(fileIds)
%CLOSEFILES Close multiple files.
%   Close any files in the list of file IDs that are actually open.

%   Copyright 1996-2003 The MathWorks, Inc.
%   $Revision: 1.1.10.3 $ $Date: 2005/11/15 01:07:12 $

for k = 1:length(fileIds)
    if ~isempty(fopen(fileIds(k)))
        fclose(fileIds(k));
    end
end
