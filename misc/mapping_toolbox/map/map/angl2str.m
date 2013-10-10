function str = angl2str(angin, signcode, units, digits)
%ANGL2STR Format angle strings
%
%   STR = ANGL2STR(ANGLE) converts a numerical vector of angles in
%   degrees to a string matrix.
%
%   STR = ANGL2STR(ANGLE, SIGNCODE) uses the string SIGNCODE to specify
%   the method for indicating that a given angle is positive or
%   negative. SIGNCODE may be one of the following:
%
%     'ew' for east/west notation
%        Trailing 'E' for positive longitudes
%        Trailing 'W' for negative longitudes
%
%     'ns' for north/south notation
%        Trailing 'N' for positive latitudes
%        Trailing 'S' for negative latitudes
%
%     'pm' for plus/minus notation
%        Leading '+' sign for positive angles
%        Leading '-' sign for negative angles
%
%     'none' for blank/minus notation
%        Sign omitted for positive angles
%        Leading '-' sign for negative angles
%
%   The default value for SIGNCODE is 'none'.
%
%   STR = ANGL2STR(ANGLE, SIGNCODE, UNITS) uses the string UNITS to
%   indicate both the units in which ANGLE is provided AND to control
%   the output format.  UNITS may be: 'degrees' (the default value),
%   'degrees2dm', 'degrees2dms', or 'radians'.  It is interpreted as
%   follows:
%
%      UNITS         Units of ANGLE    Output format
%      -----         --------------    -------------
%     'degrees'      degrees           decimal degrees
%     'degrees2dm'   degrees           degrees/decimal minutes
%     'degrees2dms'  degrees           degrees/minutes/decimal seconds
%     'radians'      radians           decimal radians        
%
%   STR = ANGL2STR(ANGLE, SIGNCODE, UNITS, N) uses the integer N to
%   control the number of significant digits provided in the output.
%   N is the power of 10 representing the last place of significance in
%   the number of degrees, minutes, seconds, or radians -- for UNITS of
%   'degrees', 'degrees2dm', 'degrees2dms', and 'radians', respectively.
%   For example, if N = -2 (the default), ANGL2STR rounds to the nearest
%   hundredth. If N = 0, ANGL2STR rounds to the nearest integer.  And if
%   N == 1, ANGL2STR rounds to the tens place, although positive values
%   of N are of little practical use.  In all cases, the interpretation
%   of the parameter N is consistent between ANGL2STR and ROUNDN.

% Copyright 1996-2007 The MathWorks, Inc.
% $Revision: 1.14.4.7 $  $Date: 2007/11/09 20:23:09 $

error(nargchk(1, 4, nargin, 'struct'))

switch(nargin)
    case 1
        signcode = 'none';
        units  = 'degrees';
        digits = -2;
    case 2
        units  = 'degrees';
        digits = -2;
    case 3
        % When we have three input arguments, the third one could be
        % either UNITS or DIGITS (N).
        if ischar(units)
            units = checkunits(units);
            digits = -2;
        else
            digits = units;
            units  = 'degrees';
        end
    case 4
        units = checkunits(units);
end

angin = ignoreComplex(angin, mfilename, 'ANGLE');

str = doAngle2str(angin, signcode, units, digits);

%-----------------------------------------------------------------------

function strout = doAngle2str(angin, signcode, units, digits)

%  Ensure that inputs are a column vector
angin = angin(:);

switch units
    case 'degrees'
        middle = formatDegrees(angin, digits);

    case 'degrees2dm'
        middle = formatDM(angin, digits);

    case 'degrees2dms'
        middle = formatDMS(angin, digits);

    case 'radians'
        middle = formatRadians(angin, digits);
end

[prefix, suffix] = buildPrefixAndSuffix(angin, signcode);

strout = buildStringMatrix(prefix, middle, suffix);

%-----------------------------------------------------------------------

function [prefix, suffix] = buildPrefixAndSuffix(angin, signcode)

%  Compute the prefix and suffix matrices.
%  Note that the * character forces a space in the output

switch lower(signcode)
   case 'ns'
      prefix = [];
      suffix = '**';
      suffix = suffix(ones(size(angin)),:);
      suffix(angin > 0, 2) = 'N';
      suffix(angin < 0, 2) = 'S';

   case 'ew'
      prefix = [];
      suffix = '**';
      suffix = suffix(ones(size(angin)),:);
      suffix(angin > 0, 2) = 'E';
      suffix(angin < 0, 2) = 'W';

   case 'pm'
      prefix = ' ';
      prefix = prefix(ones(size(angin)),:);
      prefix(angin > 0) = '+';
      prefix(angin < 0) = '-';
      suffix = [];

   case 'none'
      prefix = ' ';
      prefix = prefix(ones(size(angin)),:);
      prefix(angin < 0) = '-';
      suffix = [];

   otherwise
      eid = sprintf('map:%s:unknownFormatString', mfilename);
      error(eid,'%s','Unrecognized SIGNCODE string')

end

%-----------------------------------------------------------------------
function strout = buildStringMatrix(prefix, middle, suffix)

strout = [prefix middle suffix];

% Right justify each row of the output matrix.  This places
% all extra spaces in the leading position.  Then strip these
% lead zeros.  Left justifying and then a DEBLANK call will
% not ensure that all strings line up.  LEADBLNK only strips
% leading blanks which appear in all rows of a string matrix,
% thereby not messing up any right justification of the string matrix.
strout = strjust(strout);
strout = leadblnk(strout,' ');

% Replace the hold characters with a space
strout(strout == '*') = ' ';

% Pad matrix with a space at front and back to avoid touching the map frame
strout = [ char(32*ones(size(strout,1),1)) strout char(32*ones(size(strout,1),1))];

%-----------------------------------------------------------------------

function units = checkunits(units)

% All comparisons here are case-insensitive. Error on 'dms' or 'dm'
% because the old encodings are no longer supported.  Then check for an
% exact match to 'degrees2dms' or 'degrees2dm'.  Then check for a full
% or partial match to 'degrees' or 'radians'.

if strcmpi(units, 'dms')
    error('map:angl2str:obsoleteDMS', ...
        ['DMS angle encoding is obsolete.\n', ...
        'To format output in degrees-minutes-seconds, set UNITS to ''degrees2dms''.'])
elseif strcmpi(units, 'dm')
    error('map:angl2str:obsoleteDM', ...
        ['DM angle encoding is obsolete.\n', ...
        'To format output in degrees-minutes, set UNITS to ''degrees2dm''.'])
elseif any(strcmpi(units, {'degrees2dms','degrees2dm'}))
    units = lower(units);
else
    units = checkangleunits(units);
end

%-----------------------------------------------------------------------

function middle = formatDMS(angin, digits)

% Tony's trick
onesize = ones(size(angin));

spacestr    = '*';
spacestr    = spacestr(onesize);

quotestr    = '''';
quotestr    = quotestr(onesize);

dblquotestr = '"';
dblquotestr = dblquotestr(onesize);

%  Convert to matrix format.  Work with only the absolute values
%  since the sign is taken care of by the prefix string
dms = degrees2dms(abs(angin));

d = dms(:,1);
m = dms(:,2);
s = dms(:,3);

% Avoid the equivalent of "m = -0;" because it may cause an
% unnecessary '-' sign in the output of num2str.

m(m == 0) = abs(m(m == 0));

if digits >= 2
    % We are rounding seconds away completely;
    % make sure to go to the nearest minute.
    m = m + round(s/60);
    s = 0;
else
    s = roundn(s,digits);
end

%  Degrees, minutes and seconds

d_str = num2str(d,'%4g');        %  Convert degrees to a string
m_str = num2str(m,'%02g');       %  Convert minutes to a string
s_str = num2str(s,formatstr(digits));  %  Convert seconds to a padded string

%  Construct the display string

middle = [leadblnk(d_str) degsymbol(size(angin)) ...
    spacestr m_str quotestr spacestr s_str dblquotestr];
                    
%-----------------------------------------------------------------------

function middle = formatDM(angin, digits)

% Tony's trick
onesize = ones(size(angin));

spacestr    = '*';
spacestr    = spacestr(onesize);

quotestr    = '''';
quotestr    = quotestr(onesize);

%  Convert to matrix format.  Work with only the absolute values
%  since the sign is taken care of by the prefix string

dm = degrees2dm(abs(angin));
d = dm(:,1);
m = dm(:,2);

if digits >= 2
    % We are rounding minutes away completely;
    % make sure to go to the nearest degree.
    d = d + round(m/60);
    m = 0;
else
    m = roundn(m,digits);
end

d_str = num2str(d,'%4g');
m_str = num2str(m,formatstr(digits));

middle = [leadblnk(d_str) degsymbol(size(angin)) spacestr  m_str quotestr];

%-----------------------------------------------------------------------                   

function fmt = formatstr(digits)
%  Construct the format string for converting seconds (or minutes)

rightdigits  = abs(min(digits,0));
if rightdigits > 0
    totaldigits = 3 + rightdigits;
else
    totaldigits = 2 + rightdigits;
end

fmt = ['%0',num2str(totaldigits),'.',num2str(rightdigits),'f'];

%-----------------------------------------------------------------------

function middle = formatDegrees(angin, digits)

angin = roundn(angin,digits);
formatstr = ['%20.',num2str(abs(min(digits,0)) ),'f'];
str = num2str(abs(angin),formatstr);
middle = [leadblnk(str) degsymbol(size(angin))];
    
%-----------------------------------------------------------------------

function middle = formatRadians(angin, digits)

angin = roundn(angin,digits);
formatstr = ['%20.',num2str(abs(min(digits,0)) ),'f'];
str = num2str(abs(angin),formatstr);
unitsymbol = '*R';
unitsymbol = unitsymbol(ones(size(angin)),:);
middle = [leadblnk(str) unitsymbol];

%-----------------------------------------------------------------------

function unitsymbol = degsymbol(sizein)

unitsymbol = degchar;
unitsymbol = unitsymbol(ones(sizein),:);
