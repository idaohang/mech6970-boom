function angles = str2angle(strings)
%STR2ANGLE Convert strings to angles in degrees
%
%   ANGLES = STR2ANGLE(STRINGS) converts strings containing latitudes
%   and/or longitudes, expressed in one of four different formats with
%   degree-minutes-seconds, to numeric angles in units of degrees.
%
%      Format Description                      Example
%      ------------------                     -------
%      Degree Symbol, Single/Double Quotes    '123°30''00"W'
%      'd', 'm', 's' Separators               '123d30m00sW'
%      Minus Signs as Separators              '123-30-00W'
%      "Packed DMS"                           '1233000W'
%
%   Note: In the first example, the character that immediately follows
%         '123' is meant to be a degree symbol, but may display
%         differently depending on your locale.  If this is the case,
%         you should try to work with one of the other three formats.
%
%   Input must conform closely to the examples provided; in particular the
%   seconds field must be included, even if it is not significant.  Except
%   in Packed DMS format, the seconds field may contain a fractional
%   component.  Sign characters are not supported; terminate each entry
%   with 'N' for positive latitude, 'S' for negative latitude, 'E' for
%   positive longitude, 'W' for negative longitude.  STRINGS is a string
%   or a cell array of strings.  For backward compatibility, STRINGS can
%   also be a character matrix. If more than one angle is represented,
%   STRINGS can contain either homogeneous or heterogeneous formatting.
%   ANGLES is a column vector of class double.

% Copyright 1996-2008 The MathWorks, Inc.
% $Revision: 1.6.4.6 $  $Date: 2008/05/14 22:01:13 $

strings = preprocess(strings);

% Attempt to identify format based on the existence of certain
% characters in the string matrix. Then attempt vector processing for
% the most plausible format.  If an inconsistency is detected, process
% row-by-row instead.
if containsQuotes(strings)
   % ddd ° mm " ss.sss [N | S] , used by NIMA
   % Enhancement: Support leading [+ | -] in place of trailing
   %    hemisphere string
   % Don't actually check for the degrees symbol, because the 
   % character code varies across platforms
    angles = degreeCharacterAndQuotes(strings);
elseif containsdms(strings)
    % 111d22m33sS    
    angles = dmsSeparators(strings);
elseif usesMinusSignsAsSeparators(strings)
    % [d]dd-mm-ss[.ss](N | S), used by the FAA
    angles = minusSignsAsSeparators(strings);
elseif charSetFitsPackedDMS(strings)
    % [d]ddmmss(N | S)
    angles = packedDMS(strings);
else
    error('map:str2angle:invalidString', ...
        'Unrecognized or invalid angle format.')
    % Enhancement: Support census [+-][d]ddmmsssss
    % Enhancement: Support ERS header (like in new etopo5) [-][d]dd:mm:ss
end

%---------------------------------------------------------------------------

function angles = packedDMS(strings)

s = str2num(strings(:,end-2:end-1));   %#ok
m = str2num(strings(:,end-4:end-3));   %#ok
d = str2num(strings(:,max(1,end-7):end-5));   %#ok
if consistentColumnLength(strings, d, m, s)
    angles = hemispheresign(strings) .* (d + (m + s/60)/60);
else
    angles = processRowByRow(strings);
end   

%---------------------------------------------------------------------------

function angles = minusSignsAsSeparators(strings)

[i,j] = ind2sub(size(strings),find(strings(:)'=='-'));

inconsistentFormats = (size(strings,1) > 1) ...
    && (~all(diff(j(end/2+1:end)))==0 || ~all(diff(j(1:end/2)))==0);

if ~inconsistentFormats
    s = str2num(strings(:,j(end)+1:end-1));   %#ok
    m = str2num(strings(:,j(1)+1:j(end)-1));   %#ok
    d = str2num(strings(:,1:j(1)-1));   %#ok
    if consistentColumnLength(strings, d, m, s)
        angles = hemispheresign(strings) .* (d + (m + s/60)/60);
    else
        angles = processRowByRow(strings);
    end
else
    angles = processRowByRow(strings);
end

%---------------------------------------------------------------------------

function angles = degreeCharacterAndQuotes(strings)

[ip,jp] = ind2sub(size(strings),find(strings(:)'==''''));
[ipp,jpp] = ind2sub(size(strings),find(strings(:)'=='"'));

inconsistentFormats = ...
    (size(strings,1) > 1) && (~all(diff(jp))==0 ||  ~all(diff(jpp))==0);

if ~inconsistentFormats
    s = str2num(strings(:,jp(1)+1:jpp(1)-1));   %#ok
    m = str2num(strings(:,jp(1)-2:jp(1)-1));   %#ok
    d = str2num(strings(:,1:jp(1)-4));   %#ok
    if consistentColumnLength(strings, d, m, s)
        angles = hemispheresign(strings) .* (d + (m + s/60)/60);
    else
        angles = processRowByRow(strings);
    end
else
    angles = processRowByRow(strings);
end

%---------------------------------------------------------------------------

function angles = dmsSeparators(strings)

[id,jd] = ind2sub(size(strings),find(strings(:)'=='d'));
[im,jm] = ind2sub(size(strings),find(strings(:)'=='m'));
[is,js] = ind2sub(size(strings),find(strings(:)'=='s'));

inconsistentFormats = ...
    (size(strings,1) > 1) && (~all(diff(jd))==0 ...
    || ~all(diff(jm))==0 || ~all(diff(js))==0);

if ~inconsistentFormats
    s = str2num(strings(:,jm(1)+1:js(1)-1));   %#ok
    m = str2num(strings(:,jd(1)+1:jm(1)-1));   %#ok
    d = str2num(strings(:,1:jd(1)-1));   %#ok
    if consistentColumnLength(strings, d, m, s)
        angles = hemispheresign(strings) .* (d + (m + s/60)/60);
    else
        angles = processRowByRow(strings);
    end
else
    angles = processRowByRow(strings);
end

%---------------------------------------------------------------------------

function sgn = hemispheresign(strings)

sgn = ones(size(strings,1),1);
sgn( strings(:,end) == 'S' | strings(:,end) == 'W' ) = -1;

%---------------------------------------------------------------------------

function angles = processRowByRow(strings)
angles = zeros(size(strings,1),1);
for k = 1:size(strings,1)
    angles(k) = str2angle(strings(k,:));
end

%---------------------------------------------------------------------------

function strings = rightjustify(strings)
% Right justify a 2-D character array.
strings = fliplr(char(strtrim(cellstr(fliplr(char(strings))))));

%---------------------------------------------------------------------------

function tf = containsQuotes(strings)
tf = any(strings(:) == '''') && any(strings(:) == '"');

%---------------------------------------------------------------------------

function tf = containsdms(strings)
tf = any(strings(:) == 'd') && ...
     any(strings(:) == 'm') && ...
     any(strings(:) == 's');

%---------------------------------------------------------------------------

function tf = usesMinusSignsAsSeparators(strings)
tf = (numel(findstr(strings(:)','-')) == 2*size(strings,1));

%---------------------------------------------------------------------------

function tf = charSetFitsPackedDMS(strings)
permissibleCharacters = '0123456789NSEW+-. ';
tf = isempty(setdiff(strings(:)', double(permissibleCharacters)));

%---------------------------------------------------------------------------

function tf = consistentColumnLength(strings, d, m, s)
tf = isequal( size(strings,1), length(d), length(m), length(s) );

%---------------------------------------------------------------------------

function strings = preprocess(strings)

% convert cell array to a cell array vector
if iscell(strings)
   strings = strings(:);
end

% If the string came from the ANGL2STR, replace the Latex degree
% character substring.
if ischar(strings)
    strings = num2cell(strings,2);
end
strings = strrep(strings,'^{\circ}','°');

% Convert to character array, justify, validate
strings = str2mat(strings{:});
strings = rightjustify(strings);

assert(isempty(setdiff(strings(:,end),double('NSEW'))), ...
    'map:str2angle:illegalHemisphereCharacter', ...
    'Last character not N,S,E or W')

assert(size(strings,2) >= 7, ...
    'map:str2angle:notEnoughNumerals',...
    'Unrecognized angle format: string seems too short.')

assert(all(strings(:,1) ~= '-') && all(strings(:,1) ~= '+'), ...
            'map:str2angle:leadingSignNotAllowed', ...
            'Invalid leading sign character(s).')
