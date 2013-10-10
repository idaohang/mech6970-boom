function R = worldfileread(worldFileName, coordinateSystemType, rasterSize)
%WORLDFILEREAD Read world file and return referencing object or matrix
%
%   R = WORLDFILEREAD(worldFileName, coordinateSystemType, rasterSize)
%   reads the world file, worldFileName, and constructs a spatial
%   referencing object, R. The type of referencing object is determined by
%   the coordinateSystemType string, which can be either 'planar'
%   (including projected map coordinate systems) or 'geographic' (for
%   latitude-longitude systems). The rasterSize input should match to the
%   size of the image corresponding to the world file.
%
%   REFMAT = WORLDFILEREAD(worldFileName) reads the world file, 
%   worldFileName, and constructs a 3-by-2 referencing matrix, REFMAT.
%
%   Example 1
%   ---------
%   % Read ortho image referenced to a projected coordinate system
%   % (Massachusetts State Plane Mainland)
%   filename = 'concord_ortho_w.tif';
%   [X, cmap] = imread(filename);
%   worldFileName = getworldfilename(filename);
%   R = worldfileread(worldFileName, 'planar', size(X))
%
%   Example 2
%   ---------
%   % Read image referenced to a geographic coordinate system
%   filename = 'boston_ovr.jpg';
%   RGB = imread(filename);
%   worldFileName = getworldfilename(filename);
%   R = worldfileread(worldFileName, 'geographic', size(RGB))
%
%   See also GETWORLDFILENAME, PIX2MAP, MAP2PIX, WORLDFILEWRITE

% Copyright 1996-2010 The MathWorks, Inc.
% $Revision: 1.1.10.6 $  $Date: 2010/10/11 14:47:50 $

error(nargchk(1, 3, nargin, 'struct'))
if nargin == 2
    error('map:worldfile:expected1or3Inputs', ...
        'Function %s requires exactly 1 input or 3 inputs.', ...
        'WORLDFILEREAD')
end

if nargin == 3
    % Validate coordinateSystemType, but let georasterref or maprasterref
    % validate rasterSize.
    coordinateSystemType = validatestring( ...
        coordinateSystemType, {'geographic', 'planar'}, ...
        'WORLDFILEREAD', 'coordinateSystemType', 2);
end

% Check worldFileName must be a string
assert(ischar(worldFileName), ...
    'map:validate:expectedString', ...
    'Input argument %s must be a string.', 'WORLDFILENAME')

% Verify that the worldFileName is a file
assert(exist(worldFileName, 'file') == 2, ...
    'map:worldfile:fileDoesNotExist', ...
    'The world file with name ''%s'' does not exist.', worldFileName)

% Try to open the input worldFileName
fid = fopen(worldFileName);
if fid == -1
    error('map:worldfile:unableToOpenFile', ...
        'Unable to open file ''%s''.', worldFileName);
end
clean = onCleanup(@() fclose(fid));

% Read W into a 6-element vector
[W, count] = fscanf(fid,'%f',6);
assert(count == 6, ...
    'map:worldfile:expectedSixNumbers', ...
    'Unexpected world file contents.');

if nargin == 3
    if strcmp(coordinateSystemType,'geographic')
        R = georasterref(W, rasterSize, 'cells');
    else
        R = maprasterref(W, rasterSize, 'cells');
    end
else
    % Convert W to a referencing matrix
    R = worldFileMatrixToRefmat(W);
end
