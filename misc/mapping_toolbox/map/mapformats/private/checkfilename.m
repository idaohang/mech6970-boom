function varargout = checkfilename(varargin)
%CHECKFILENAME Check validity of a filename.
%   FULLFILENAME = CHECKFILENAME( ...
%       FILENAME, FUNCTION_NAME, ARGUMENT_POSITION)
%   checks the validity of the FILENAME for reading data and issues a
%   formatted error if the FILENAME is invalid. FULLFILENAME is the
%   absolute pathname of the file.
%
%   FILENAME can be a MATLABPATH relative partial pathname. If the file is
%   not found in the current working directory, CHECKFILENAME searches down
%   MATLAB's search path. The FILENAME must exist with read permission. If
%   FILENAME is a URL address, an error is issued.
%
%   FUNCTION_NAME is a string containing the function name to be used in
%   the formatted error message.
%
%   ARGUMENT_POSITION is a positive integer indicating which input argument
%   is being checked; it is also used in the formatted error message.
%
%   FULLFILENAME = CHECKFILENAME( ...
%       FILENAME, EXT, FUNCTION_NAME, ARGUMENT_POSITION)
%   loops, appending each element of EXT in both upper and lower case,
%   until a valid filename is found; otherwise it will issue an error.  EXT
%   is a cell array of strings, which may be EMPTY, containing extensions
%   without the '.' character.
%
%   [FULLFILENAME, FID] = CHECKFILENAME(FILENAME, ...) leaves the file
%   FULLFILENAME opened and returns the file identifier.
%
%   [FULLFILENAME, ISURL] = CHECKFILENAME(FILENAME, ..., ALLOWURL) allows
%   FILENAME to be a URL address, if ALLOWURL is true. The address must
%   include the protocol type (e.g., "http://").  In this case,
%   FULLFILENAME is the absolute path to a temporary file copied from the
%   URL location and ISURL is TRUE; otherwise, ISURL is false.  If FILENAME
%   is a URL address, remember to delete FULLFILENAME as it is a temporary
%   copy.

% Copyright 1996-2010 The MathWorks, Inc.
% $Revision: 1.1.6.4 $  $Date: 2010/03/04 16:22:29 $

error(nargoutchk(0, 2, nargout, 'struct'))

% Validate inputs and outputs.
[filename, ext, function_name, allowURL] = parseInputs(varargin{:});
haveAllowURL = (nargin >= 4) && islogical(varargin{end});

if (strfind(filename, '://'))
    % The filename is a URL.
    fullfilename=checkURL(filename, function_name, allowURL, nargout);
    fid = [];
    isURL = true;
else
    % The filename is a disk filename.
    % fid points to an open file.
    [fullfilename, fid]=checkDiskfilename(filename, ext, function_name);
    isURL = false;

    % Close the file if not requesting a FID output or
    % if using the ALLOWURL syntax.
    needToCloseFile = nargout ~= 2 || haveAllowURL;
    if needToCloseFile
        fclose(fid);
    end
end

needFID = (nargout == 2) && ~haveAllowURL;
needURL = (nargout == 2) &&  haveAllowURL;
optionalOutputs = {fid, isURL};
optionalOutputs([~needFID, ~needURL]) = [];
varargout = {fullfilename optionalOutputs{:}};

%-------------------------------------------------------------------------
function [filename, ext, function_name, allowURL] = parseInputs(varargin)

allowURL= false;
switch nargin
   case {0,1,2}
      % CHECKFILENAME
      % CHECKFILENAME(FILENAME)
      % CHECKFILENAME(FILENAME, FUNCTION_NAME)
      eid = sprintf('%s:%s:tooFewInputs', getcomp, mfilename);
      error(eid, '%s', 'Error using CHECKFILENAME: too few inputs.');

   case 3
      % CHECKFILENAME(FILENAME, FUNCTION_NAME, ARGUMENT_POSITION)
      filename = varargin{1};
      function_name = varargin{2};
      pos = varargin{3};
      ext = {''};

   case 4
      filename = varargin{1};
      if islogical(varargin{4})
         % CHECKFILENAME(FILENAME, FUNCTION_NAME, ...
         %               ARGUMENT_POSITION, ALLOWURL)
         function_name = varargin{2};
         pos = varargin{3};
         allowURL = varargin{4};
         ext = {''};
      else
         % CHECKFILENAME(FILENAME, EXT, FUNCTION_NAME, ...)
         %               ARGUMENT_POSITION)
         ext = varargin{2};
         function_name = varargin{3};
         pos = varargin{4};
      end

   case 5
      % CHECKFILENAME(FILENAME, EXT, FUNCTION_NAME, ,...
      %               ARGUMENT_POSITION, ALLOWURL)
      filename = varargin{1};
      if islogical(varargin{5})
         ext = varargin{2};
         function_name = varargin{3};
         pos = varargin{4};
         allowURL = varargin{5};
      else
         eid = sprintf('%s:%s:invalidInputs', getcomp, mfilename);
         error(eid, '%s', ...
            'Error using CHECKFILENAME: invalid ALLOWURL argument.');
      end

   otherwise
      eid = sprintf('%s:%s:tooManyInputs', getcomp, mfilename);
      error(eid, '%s', ...
         'Error using CHECKFILENAME: too many input arguments.');
end

% Verify the inputs
if ~ischar(function_name)
   eid = sprintf('%s:%s:invalidFunctionName', getcomp, mfilename);
   error(eid, '%s', ...
      'Error using CHECKFILENAME: invalid FUNCTION_NAME argument.');
end
checkinput(filename, {'char'}, {'nonempty'}, function_name, 'FILENAME', pos);

if ~iscell(ext)
   eid = sprintf('%s:%s:invalidExt', getcomp, mfilename);
   error(eid, '%s', 'Error using CHECKFILENAME: invalid EXT argument.');
end

if ~isnumeric(pos)
   eid = sprintf('%s:%s:invalidPosition', getcomp, mfilename);
   error(eid, '%s', 'Error using CHECKFILENAME: invalid ARGUMENT_POSITION argument.');
end

%-------------------------------------------------------------------------
function [fullfilename, fid] = ...
   checkDiskfilename(filename, ext, function_name)

% Try to open the filename with read permission.
fid = fopen(filename, 'r');

% Verify the file is opened.
if (fid == -1)
   % The file cannot be opened with read permission.
   % Try to append various extensions and try again.
   found = 0;
   extension = [lower(ext) upper(ext)];
   for p = 1:numel(extension)
      extname = [filename '.' extension{p}];
      fid = fopen(extname, 'r');
      if (fid ~= -1)
         % File was found.
         %  Break "for" loop with updated extname.
         found = 1;
         break;
      end
   end

   % Check that some filename+extension combination was found.
   if (~found)
      % Filename could not be opened with 'r' permission
      %  but could exist with write mode or as a MATLAB file
      if isequal(exist(filename,'file'),2)
         % Check to verify NOT a MATLAB file
         if ~(isequal(exist([filename '.m'],'file'),2) || ...
               isequal(exist([filename '.M'],'file'),2))
            % File exists without read mode
            eid = sprintf('%s:%s:invalidFileMode', getcomp, function_name);
            error(eid, ...
                'Function %s was unable to open file ''%s'' for reading.', ...
                upper(function_name), filename);
         end
      end
      % File does not exist.
      eid = sprintf('%s:%s:invalidFilename', getcomp, function_name);
      error(eid,'Function %s was unable to find file ''%s''.', ...
         upper(function_name), filename);
   end
end

% Return the full pathname if not in pwd.
fullfilename = fopen(fid);
if isequal(exist(fullfile(pwd,fullfilename),'file'),2)
   fullfilename = fullfile(pwd, fullfilename);
end

%-------------------------------------------------------------------------
function fullfilename =  ...
   checkURL(filename, function_name, allowURL, fcnNargout)

if ~usejava('jvm')
   eid = sprintf('%s:%s:invalidJVM', getcomp, function_name);
   error(eid, ...
       'Function %s requires a Java Virtual Machine to read a URL.', ...
       upper(function_name));
end

if ~allowURL
   eid = sprintf('%s:%s:invalidFilename', getcomp, function_name);
   error(eid, 'File ''%s'' is a URL address.', filename);
end

if fcnNargout < 2
   error(nargoutchk(2, 2, fcnNargout, 'struct'))
end

urltempname = tempname;
try
   fullfilename = urlwrite(filename, urltempname);
catch %#ok<CTCH>
   % The PC will not create urltempname,
   % but Unix will, so delete it if it exists.
   if exist(urltempname, 'file')
      delete(urltempname);
   end
   eid = sprintf('%s:%s:invalidURL', getcomp, function_name);
   error(eid,'Function %s was unable to read URL ''%s''.', ...
      upper(function_name), filename);
end
