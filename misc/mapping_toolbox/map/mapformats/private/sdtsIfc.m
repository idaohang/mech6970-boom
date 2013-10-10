function [info, data] = sdtsIfc(filename)
% SDTSIfc Interface to the SDTS++ library.
%
%   INFO = SDTSIFC(FILENAME) returns a structure whose fields contain
%   information about the contents of an SDTS data set.
%
%   [INFO, DATA] = SDTSIFC(FILENAME) returns the INFO structure and reads
%   data from an SDTS raster or DEM data set.  DATA is a matrix containing
%   the elevation/data values.   
%
%   FILENAME is a string that specifies the name of the SDTS catalog
%   directory file, such as 7783CATD.DDF.  The FILENAME may also include
%   the directory name.  If FILENAME does not include the directory, then
%   it must be in the current directory or in a directory on the MATLAB
%   path.  
%
%   SDTSIFC is a wrapper function with no argument checking around the
%   SDTSMEX MEX-function. On Windows, the function must CD to the data
%   directory.
%    
%   Example
%   -------
%   info = sdtsIfc('9129CATD.DDF');
%
%   See also SDTSDEMREAD, SDTSINFO.

% Copyright 2005-2007 The MathWorks, Inc.
% $Revision: 1.1.6.3 $  $Date: 2007/12/10 21:39:41 $   

try
   if ispc
      cwd = pwd;
      pathstr = fileparts(filename);
      cd(pathstr)
   end
   switch nargout
      case {0,1}
         info = sdtsmex(filename);
      case 2
         [info, data]= sdtsmex(filename);
      otherwise
        error('map:sdtsIfc:internalError','Incorrect number of outputs.');
   end
   if ispc
      cd(cwd)
   end
catch e
   if ispc
      cd(cwd)
   end
   rethrow(e)
end
