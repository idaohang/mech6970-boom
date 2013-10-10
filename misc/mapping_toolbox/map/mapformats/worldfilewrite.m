function worldfilewrite(R, worldfilename)
%WORLDFILEWRITE Write world file from referencing object or matrix
%
%   WORLDFILEWRITE(R, WORLDFILENAME) calculates the world file entries
%   corresponding to referencing object or matrix R and writes them into
%   the file WORLDFILENAME. R can be a spatialref.MapRasterReference
%   object, a spatialref.GeoRasterReference object, or a 3-by-2
%   referencing matrix.
%
%   Example
%   -------
%   % Write out the information from a referencing object for the
%   % image file concord_ortho_w.tif
%   info = imfinfo('concord_ortho_w.tif');
%   R = worldfileread('concord_ortho_w.tfw', ...
%           'planar', [info.Height info.Width])
%   worldfilewrite(R, 'concord_ortho_w_test.tfw');
%   type concord_ortho_w_test.tfw
%
%   See also GETWORLDFILENAME, PIX2MAP, MAP2PIX, WORLDFILEREAD.

% Copyright 1996-2010 The MathWorks, Inc.
% $Revision: 1.1.10.7.2.1 $  $Date: 2010/12/03 21:43:53 $

% Validate input arguments
checkfilenameforwriting(worldfilename,'WORLDFILENAME')
classes = {'spatialref.MapRasterReference','spatialref.GeoRasterReference'};
internal.map.validateRasterReference(R, classes, 'WORLDFILEWRITE', 'R', 1)

% Try to open the output worldfilename
fid = fopen(worldfilename, 'w');

assert(fid ~= -1, ...
    'map:worldfilewrite:failedToOpen', ...
    'Unable to open the file ''%s'' for writing.', worldfilename);

clean = onCleanup(@() fclose(fid));

% Obtain a "world file matrix" W
if isobject(R)
    W = R.worldFileMatrix();
else
    W = refmatToWorldFileMatrix(R);
end

% Write W to the world file
for k = 1:numel(W)
    fprintf(fid, '%s\n', fullPrecisionFixedPointString(W(k)));
end

%--------------------------------------------------------------------------

function checkfilenameforwriting(filename, var_name)

% Filename must be a string
assert(ischar(filename), ...
    'map:worldfilewrite:expectedString', ...
    '%s must be a string.', var_name)

% Filename must not be not a directory
assert(~exist(filename, 'dir'), ...
    'map:worldfilewrite:filenameIsDirectory', ...
    'The input %s ''%s'' is a directory.', var_name, filename)

%--------------------------------------------------------------------------

function str = fullPrecisionFixedPointString(x)
%fullPrecisionFixedPointString Number to full precision fixed point string
%
%   Convert a real floating-point scalar X (class double) to a
%   fixed-point decimal string while retaining full precision. Remove
%   trailing zeros while leaving at least one digit to the right of the
%   decimal place.

% Try to print a string with enough decimal places such that the
% original value can be recovered using str2double.
nDigits = floor(-log10(eps(x)));
str = sprintf('%.*f', nDigits, x);

if ~isequal(str2double(str), x)
    % Add more precision only if necessary.
    str = sprintf('%.*f', nDigits + 1, x); 
end

% Remove trailing zeros while leaving at least one digit to the right of
% the decimal place.
n = regexp(str,'\.0*$');
if ~isempty(n)
    % There is nothing but zeros to the right of the decimal place;
    % the value in n is the index of the decimal place itself.
    % Remove all trailing zeros except for the first one.
    str(n+2:end) = [];
else
    % There is a non-zero digit to the right of the decimal place.
    m = regexp(str,'0*$');
    if ~isempty(m)
        % There are trailing zeros, and the value in m is the index of
        % the first trailing zero. Remove them all.
        str(m:end) = [];
    end
end
