function strout = dist2str(distin,format,units,digits)
%DIST2STR  Format distance strings
%
%  str = DIST2STR(dist) converts a numerical vector of distances in
%  kilometers to a string matrix.  The output string matrix is useful for
%  the display of distances.
%
%  str = DIST2STR(dist,'format') uses the specified format input
%  to construct the string matrix.  Allowable format strings are
%  'pm' for plus/minus notation; and 'none' for blank/minus notation.
%  If omitted or blank, 'none' is assumed.
%
%  str = DIST2STR(dist,'format','units') defines the units in which the
%  input distances are supplied, and which are encoded in the string
%  matrix.  Units must be one of the following: 'feet', 'kilometers',
%  'meters', 'nauticalmiles', 'statutemiles', 'degrees', or 'radians'.
%  Note that statute miles are encoded as 'mi' in the string matrix,
%  whereas in most Mapping Toolbox functions, 'mi' indicates international
%  miles. If omitted or blank, 'kilometers' is assumed.
%
%  str = DIST2STR(dist,'format',digits) uses the input digits to
%  determine the number of decimal digits in the output matrix.
%  n = -2 uses accuracy in the hundredths position, n = 0 uses
%  accuracy in the units position.  Default is n = -2.  For further
%  discussion of the input n, see ROUNDN.
%
%  str = DIST2STR(dist,'format','units',digits) uses all the inputs
%  to construct the output string matrix.

% Copyright 1996-2007 The MathWorks, Inc.
% $Revision: 1.10.4.5 $  $Date: 2007/11/09 20:23:20 $
% Written by:  E. Byrns, E. Brown

error(nargchk(1, 4, nargin, 'struct'))

if nargin == 1
    format = 'none';
    units  = 'kilometers';
    digits = -2;
elseif nargin == 2
    units  = 'kilometers';
    digits = -2;
elseif nargin == 3
    % The third argument could be either format or digits.
    if ischar(units)
        units  = unitstrd(units);
        digits = -2;
    else
        digits = units;
        units = 'kilometers';
    end
elseif nargin == 4
    units = unitstrd(units);
end

assert(ischar(format), ...
    ['map:' mfilename ':nonCharFormatString'], ...
    'Input FORMAT must be a string.')

%  Prevent complex distances
distin = ignoreComplex(distin,mfilename,'dist');

%  Ensure that inputs are a column vector

distin = distin(:);

%  Compute the prefix and suffix matrices.
%  Note that the * character forces a space in the output

switch lower(format)
   case 'pm'
      prefix = ' ';     prefix = prefix(ones(size(distin)),:);
      indx = find(distin>0);  if ~isempty(indx);   prefix(indx) = '+';   end
      indx = find(distin<0);  if ~isempty(indx);   prefix(indx) = '-';   end

   case 'none'
      prefix = ' ';     prefix = prefix(ones(size(distin)),:);
      indx = find(distin<0);  if ~isempty(indx);   prefix(indx) = '-';   end

   otherwise
      error(['map:' mfilename ':mapError'], 'Unrecognized format string')

end


%  Compute the units suffix

switch units
	case 'degrees',         suffix = degchar;
    case 'kilometers',      suffix = '*km';
    case 'nauticalmiles',   suffix = '*nm';
	case 'statutemiles',    suffix = '*mi';
	case 'radians',         suffix = '*R';
	case 'meters',          suffix = '*m';
	case 'feet',            suffix = '*ft';
end

%  Expand the suffix matrix to the same length as the input vector

suffix = suffix(ones(size(distin)),:);

%  Convert the distance vector to a string format

formatstr = ['%20.',num2str(abs(min(digits,0)) ),'f'];
str = num2str(abs(distin),formatstr);      %  Convert to a padded string
strout = [prefix str suffix];              %  Construct output string

%  Right justify each row of the output matrix.  This places
%  all extra spaces in the leading position.  Then strip these
%  lead zeros.  Left justifying and then a DEBLANK call will
%  not ensure that all strings line up.  LEADBLNK only strips
%  leading blanks which appear in all rows of a string matrix,
%  thereby not messing up any right justification of the string matrix.

strout = shiftspc(strout);
strout = leadblnk(strout,' ');

%  Replace the hold characters with a space

indx = find(strout == '*');
if ~isempty(indx)
    strout(indx) = ' ';
end
