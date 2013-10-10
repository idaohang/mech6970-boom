function [Z,R] = arcgridread( filename )
%ARCGRIDREAD Read gridded data set in Arc ASCII Grid Format
%
%   [Z, R] = ARCGRIDREAD(FILENAME) reads a grid from a file in Arc ASCII
%   Grid format.  Z is a 2D array containing the data values.  R is a
%   referencing matrix (see MAKEREFMAT).  NaN is assigned to elements of Z
%   corresponding to null data values in the grid file.
%
%   Example
%   -------
%   % Load and view Mount Washington terrain elevation
%   [Z,R] = arcgridread('MtWashington-ft.grd');
%   mapshow(Z,R,'DisplayType','surface');
%   xlabel('x (easting in meters)'); ylabel('y (northing in meters)')
%   demcmap(Z)
%
%   % View the terrain in 3D
%   axis normal; view(3); axis equal; grid on
%   zlabel('elevation in feet')
%
%   See also MAKEREFMAT, MAPSHOW, SDTSDEMREAD.

% Copyright 1996-2011 The MathWorks, Inc.
% $Revision: 1.1.10.6.14.1 $  $Date: 2011/01/29 14:47:39 $

% Verify the filename and check if it is a URL.
[filename, isURL] = checkfilename(filename, {'grd'}, mfilename, 1, true);

% Open the file.
fid = fopen(filename,'r');
if fid == -1
    error('map:arcgridread:internalError', ...
        'ARCGRIDREAD failed to open file ''%s''.', filename);
end

% Read the 6-line header.
hdr = readHeader(fid, isURL);

% Read the matrix of data values, putting the k-th row in the data
% file into the k-th column of matrix Z.  Close file -- nothing left to
% read after this.
Z = fscanf(fid,'%g',[hdr.ncols,hdr.nrows]);

% Close the fileID.
closeFileID(fid, isURL);

% Replace each no-data value with NaN.
Z(Z == hdr.nodata_value) = NaN;

% Orient the data so that rows are parallel to the x-axis and columns
% are parallel to the y-axis (for compatibility with MATLAB functions
% like SURF and MESH).
Z = Z';

% Construct the referencing matrix.
R = makerefmat(hdr.xllcorner + hdr.cellsize/2,...
    hdr.yllcorner + (hdr.nrows - 1/2) * hdr.cellsize,...
    hdr.cellsize, -hdr.cellsize);

%--------------------------------------------------------------------------

function hdr = readHeader(fid, isURL)

% The file header contains the following meta-data tags:
hdr = struct( ...
    'ncols','', ...
    'nrows','', ...
    'xllcorner','', ...
    'yllcorner','',...
    'cellsize','', ...
    'nodata_value','');

% Read the header line-by-line.
tags = fieldnames(hdr);
for k = 1:numel(tags)
    line = fgetl(fid);
    [tag, value] = strtok(line);
    tag = lower(tag);
    
    % Verify the tag matches an expected header tag.
    if ~any(strcmpi(tag,tags))
        closeFileID(fid, isURL);
        error('map:arcgridread:unexpectedItemInHeader', ...
            'Unexpected tag ''%s'' in file header (line %d).',...
            tag, k);
    end  
    hdr.(tag) = str2double(value);
end

%--------------------------------------------------------------------------

function closeFileID(fid, isURL)

% Obtain the filename from the file ID.
filename = fopen(fid);

% Close the file.
fclose(fid);

% If the file was downloaded from a URL, then delete the temporary file.
if (isURL)
    deleteDownload(filename);
end
