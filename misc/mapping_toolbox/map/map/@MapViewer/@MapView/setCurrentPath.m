function setCurrentPath(this,filename)

%   Copyright 1996-2003 The MathWorks, Inc.
%   $Revision: 1.1.6.3 $  $Date: 2009/12/31 18:50:48 $

pathstr = fileparts(filename);
if isempty(pathstr)
  f = which(filename);
  pathstr = fileparts(f);
end
this.CurrentPath = [pathstr filesep];

