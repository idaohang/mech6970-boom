function S = gshhs(varargin)
%GSHHS Read Global Self-consistent Hierarchical High-resolution Shoreline
%
%   S = GSHHS(FILENAME) reads GSHHS vector data for the entire world from
%   FILENAME.  GSHHS files must have names of the form 'gshhs_x.b',
%   'wdb_borders_x.b', or 'wdb_rivers_x.b' where x is one of the letters c,
%   l, i, h  and f, corresponding to increasing resolution (and file size).
%   The result returned in S is a polygon or line geographic data structure
%   array (a "geostruct", with 'Lat' and 'Lon' coordinate fields).
%
%   S = GSHHS(FILENAME, LATLIM, LONLIM) subsets the data from FILENAME
%   using the quadrangle defined by the latitude and longitude limits in
%   LATLIM and LONLIM. Any feature whose bounding box intersects the
%   quadrangle is returned, without trimming.  LATLIM and LONLIM are
%   2-element vectors with latitude and longitude in degrees in ascending
%   order. Longitude limits range from [-180 195]. If LATLIM is empty the
%   latitude limits are [-90 90]. If LONLIM is empty, the longitude limits
%   are [-180 195].
%
%   INDEXFILENAME = GSHHS(FILENAME,'createindex') creates an index file for
%   faster data access when requesting a subset of a larger dataset. The
%   index file has the same name as the GSHHS data file, but with the
%   extension 'i', instead of 'b' and is written in the same directory as
%   FILENAME. The name of the index file is returned, but no coastline data
%   are read. A call using this option should be followed by an additional
%   call to GSHHS to import actual data.
%
%   Output structure
%   ----------------
%   The output structure S contains the following fields.
%   All latitude and longitude values are in degrees.
%
%     Field name            Field contents
%     ----------            -----------------
%     'Geometry'            'Polygon' or 'Line'
%
%     'BoundingBox'         [min Lon min Lat;
%                            max Lon max Lat]
%
%     'Lon'                 Coordinate vector
%
%     'Lat'                 Coordinate vector
%
%     'South'               Southern latitude boundary
%
%     'North'               Northern latitude boundary
%
%     'West'                Western longitude boundary
%
%     'East'                Eastern longitude boundary
%
%     'Area'                Area of polygon in square kilometers
%
%     'Level'               Scalar value ranging from 1 to 4,
%                           indicates level in topological hierarchy
%
%     'LevelString'         'land', 'lake', 'island_in_lake',
%                           'pond_in_island_in_lake', or ''
%
%     'NumPoints'           Number of points in the polygon
%
%     'FormatVersion'       Format version of data file
%                           Positive integer for versions 3 and later,
%                           empty for versions 1 and 2
%
%     'Source'              Source of data: 'WDBII' or 'WVS'
%
%     'CrossGreenwich'      Scalar flag: true if the polygon crosses
%                           the prime meridian, false otherwise
%
%     'GSHHS_ID'            Unique polygon scalar ID number, starting at 0
%
%   For releases 2.0 and higher (FormatVersion 7 and higher) the following
%   additional fields are included in the output structure:
%
%     'RiverLake'           Scalar flag: true if the polygon is the fat
%                           part of a major river and level value set to 2;
%                           false otherwise.
%
%     'AreaFull'            Area of original full-resolution polygon in
%                           units 1/10 km^2
%
%     'Container'           ID of container polygon that encloses this
%                           polygon and set to -1 to indicate none
%
%     'Ancestor'            ID of ancestor polygon in the full resolution
%                           set that was the source of this polygon and
%                           set to -1 to indicate none
%
%   WDB Data Files
%   --------------
%   You can read the WDB rivers and borders datasets but the LevelString
%   field will be empty. The Level values vary from feature to feature but
%   the intepretations of these values are not documented as part of the
%   GSHHS distribution and are therefore not converted to strings.
%
%   Source
%   ------
%   For details on locating GSHHS data for download over the Internet, see
%   the following documentation at the MathWorks web site:
%
%   <a href="matlab:
%   web('http://www.mathworks.com/support/tech-notes/2100/2101.html#gshhs')
%   ">http://www.mathworks.com/support/tech-notes/2100/2101.html</a>
%
%   Example 1
%   ---------
%   % Read the entire coarse data set
%   % and display as a coastline.
%   filename = gunzip('gshhs_c.b.gz', tempdir);
%   world = gshhs(filename{1});
%   delete(filename{1})
%   figure
%   worldmap world
%   geoshow([world.Lat], [world.Lon])
%
%   % Display each level in a different color.
%   levels = [world.Level];
%   land = (levels == 1);
%   lake = (levels == 2);
%   island = (levels == 3);
%   figure
%   worldmap world
%   geoshow(world(land),  'FaceColor', [0 1 0])
%   geoshow(world(lake),  'FaceColor', [0 0 1])
%   geoshow(world(island),'FaceColor', [1 1 0])
%
%   Example 2
%   ---------
%   % Read and display Africa as a green polygon.
%   filename = gunzip('gshhs_c.b.gz', tempdir);
%   indexname = gshhs(filename{1}, 'createindex');
%   figure
%   worldmap Africa
%   projection = gcm;
%   latlim = projection.maplatlimit;
%   lonlim = projection.maplonlimit;
%   africa = gshhs(filename{1}, latlim, lonlim);
%   delete(filename{1})
%   delete(indexname)
%   geoshow(africa, 'FaceColor', 'green')
%   setm(gca, 'FFaceColor', 'cyan')
%
%   See also GEOSHOW, SHAPEREAD, WORLDMAP.

% Copyright 1996-2011 The MathWorks, Inc.
% $Revision: 1.1.6.9.8.1 $  $Date: 2011/02/22 03:33:45 $

% Check number of inputs.
checknargin(1,3,nargin,mfilename);

% Parse the inputs.
[latlim, lonlim, filename, createindex, subset, isWDB] = parseInputs(varargin);

% Create the index file, if requested.
if createindex
    ifilename = createIndexFile(filename, isWDB);
    S = strrep(ifilename,'\','/');
else
    % Get the index filename if it exists.
    ifilename = getIndexFilename(filename);
    
    % Read the data file.
    if ~isempty(ifilename) && subset
        % An index file exists. Read only the parts within the limits.
        S = gshhsReadUsingIndex(filename, ifilename, isWDB, latlim, lonlim);
    else
        % Read all of the file keeping only parts within coordinate limits.
        S = gshhsRead(filename, isWDB, latlim, lonlim, subset);
    end
end

%--------------------------------------------------------------------------

function [latlim, lonlim, filename, createindex, subset, isWDB] = parseInputs(inputs)
% Parse and return input arguments.

dataLatLim = [-90 90];
dataLonLim = [-180 195];
createindex = false;

% Input checks.
switch numel(inputs)
    case 1
        % gshhs(filename)
        filename = inputs{1};
        subset = 0;
        latlim = dataLatLim;
        lonlim = dataLonLim;
        
    case 2
        % gshhs(filename, 'createindex')
        filename = inputs{1};
        param = 'createindex';
        validatestring(inputs{2}, {param}, upper(mfilename), param, 2);
        createindex = true;
        subset = 0;
        latlim = dataLatLim;
        lonlim = dataLonLim;
        
    case 3
        % gshhs(filename, latlim, lonlim)
        filename = inputs{1};
        latlim = inputs{2};
        lonlim = inputs{3};
        if isempty(latlim)
            latlim = dataLatLim;
        end
        if isempty(lonlim)
            lonlim = dataLatLim;
        end
        subset = 1;
end

assert(ischar(filename), ...
    'map:gshhs:invalidType', 'FILENAME must be a char input.');

assert(exist(filename,'file') == 2, ...
    'map:gshhs:fileNotFound', ...
    'GSHHS data file "%s" does not exist.',filename);

[~, fname] = fileparts(filename);
isGSHHS = ~isempty(strfind(fname, 'gshhs'));
isWDB = ~isempty(strfind(fname, 'rivers')) || ~isempty(strfind(fname, 'borders'));
assert(isGSHHS || isWDB, ...
    'map:gshhs:incorrectFilenameForm', ...
    'The GSHHS data filename "%s" must be of the form: %s', ...
    filename, 'gshhs_x.b, wdb_rivers_x.b, or wdb_borders_x.b');
 
checkinput(latlim,{'double'},{'real','vector','finite'},mfilename, ...
    'LATLIM',2);
checkinput(lonlim,{'double'},{'real','vector','finite'},mfilename, ...
    'LONLIM',3);

assert(numel(latlim) == 2, ...
    'map:gshhs:invalidLatLim', ...
    'LATLIM must be a two element vector in units of degrees.');

assert(numel(lonlim) == 2, ...
    'map:gshhs:invalidLonLim', ...
    'LONLIM must be a two element vector in units of degrees.');

assert(latlim(1) < latlim(2), ...
    'map:gshhs:invalidLatOrder', ...
    '1st element of LATLIM must be greater than 2nd.')

assert(lonlim(1) < lonlim(2), ....
    'map:gshhs:invalidLonOrder', ...
    '1st element of LONLIM must be greater than 2nd.')

assert(-90 <= latlim(1) && latlim(2) <= 90, ...
    'map:gshhs:invalidLatLimits', ...
    'LATLIM must be between -90 and 90.')

assert(all(dataLonLim(1) <= lonlim) && all(lonlim <= dataLonLim(2)), ...
    'map:gshhs:invalidLonLimits', ...
    'LONLIM must be between %g and %g.', dataLonLim(1), dataLonLim(2))

%--------------------------------------------------------------------------

function ifilename = createIndexFile(filename, isWDB)
% Create the GSHHS index file.

% Open the GSHHS data file.
FileID = fileopen(filename);

% Open the index file.
ifilename = filename;
ifilename(end) = 'i';
iFileID = fopen(ifilename,'w','ieee-be');
if iFileID == -1
    fclose(FileID);
    error('map:gshhs:invalidIndexFileModes', ...
        'Unable to open index file ''%s'' for writing.', ifilename)
end

% Get the end and beginning file positions.
[EOF, BOF] = getEofBofFilePositions(FileID);

% For each polygon, read header block, and if within limits read the
% coordinates.
degreesPerMicroDegree = 1.0E-06;
dataBlockLengthInBytes = 8; % Lat, Lon in INT32
FilePosition = BOF;

% Extract the version number.
version = extractVersion(FileID);

% Obtain a handle to a reader function based on the version number.
[readHeaderFcn, extraAttributes] = getReaderFcn(FileID, version);

% Pre-allocate space for the geostruct.
S = initGeoStruct(1, extraAttributes);

% Read the data from the GSHHS file and write the number of coordinate
% points bounding box to the index file.
while FilePosition ~= EOF  
    % Read header info.
    S = readAttributes(FileID, readHeaderFcn, S, isWDB);
    
    % Write this record to the index file.
    record = [S.NumPoints [S.West S.East S.South S.North]./degreesPerMicroDegree];
    fwrite(iFileID, record,'int32');
    
    % Move to the end of this data block.
    offset = S.NumPoints * dataBlockLengthInBytes;
    FilePosition = fileseek(FileID, offset, 'cof');   
end

% Close the files,
fclose(FileID);
fclose(iFileID);

%--------------------------------------------------------------------------

function ifilename = getIndexFilename(filename)
% Return the name of the index file.

ifilename = filename;
ifilename(end) = 'i';
if ~exist(ifilename,'file')
    ifilename(end) = 'I';
    if ~exist(ifilename,'file')
        ifilename = [];
    end
end

%--------------------------------------------------------------------------

function S = gshhsReadUsingIndex(filename, ifilename, isWDB, latlim, lonlim)
% Read the GSHHS file using an index file to improve performance.

% Open the GSHHS file.
FileID = fileopen(filename);

% Open the index file again and get the number of records.
iFileID = fopen(ifilename,'rb','ieee-be');
if iFileID == -1
    fclose(FileID);
    error('map:gshhs:indexFopenError', ...
        'Unable to open index file ''%s'' for reading.',ifilename)
end
fileseek(iFileID, 0, 'eof');
ifilelength = ftell(iFileID);
fclose(iFileID);

% Use the index file and determine the file record indices of the polygons
% that are within the coordinate limits.
[extractindx, npolypts] = inlimitpolys(ifilename, ifilelength, latlim, lonlim);

if ~isempty(extractindx)     
    % Extract the version number
    version = extractVersion(FileID);
    
    % Read the data using a reader based on the version number.
    [readHeaderFcn, extraAttributes] = getReaderFcn(FileID, version);
    
    % Initialize the geostruct, S to size [length(extractindex),1].
    S = initGeoStruct(length(extractindx), extraAttributes);
    
    % Get the number of header bytes.
    numberOfHeaderBytes = getNumberOfHeaderBytes( ...
        FileID, readHeaderFcn, S(1), isWDB);
    
    % For each polygon within limits, read header block and coordinates.
    for k=1:numel(extractindx)
        % Read header info. Skip data from preceding polygons.
        bytesbefore = (extractindx(k)-1)*numberOfHeaderBytes + ...
            sum( npolypts(1:extractindx(k)-1,1) )*8;
        if ~isempty(bytesbefore)
            fileseek(FileID, bytesbefore, 'bof');
        end
        
        % Read the attribute header.
        S(k) = readAttributes(FileID, readHeaderFcn, S(k), isWDB);
        
        % Read the coordinate data.
        [S(k).Lat, S(k).Lon, S(k).BoundingBox] = readData(FileID, S(k));
    end
else
    % There are no polygons within the specified limits.
    % Initialize an empty geostruct.
    S = initGeoStruct(1);
    S(1) = [];
end

% Close the file.
fclose(FileID);

%--------------------------------------------------------------------------

function S = gshhsRead(filename, isWDB, latlim, lonlim, subset)
% Read records in a GSHHS file and return those that are within the
% coordinate limits.

% Open the GSHHS file.
FileID = fileopen(filename);

% Get the end and beginning file positions and validate the file contains
% data.
EOF = getEofBofFilePositions(FileID);

% Extract the version number.
version = extractVersion(FileID);

% For each polygon, read header block,
% and if within limits read the coordinates.
dataBlockLengthInBytes = 8;

% Read the data using a reader based on the version number.
[readHeaderFcn, extraAttributes] = getReaderFcn(FileID, version);

% Make sure the file pointer is at the beginning of the file.
fileseek(FileID, 0, 'bof');
FilePosition = ftell(FileID);

% For performance, pre-allocate a fixed number of elements of S.
blockSize = 1000;
S = initGeoStruct(blockSize, extraAttributes);
k = 0;

% Read all records until the end of the file.
while FilePosition ~= EOF
    k = k + 1;
    if k > numel(S)
        % Add blockSize more elements to S.
        S = [S; initGeoStruct(blockSize, extraAttributes)];%#ok<AGROW>
    end
    
    % Read header info.
    S(k) = readAttributes(FileID, readHeaderFcn, S(k), isWDB);
    
    % Determine if any of the data in the current block falls within the
    % input limits. If so, read the coordinates.  If not, move to the end
    % of the record.
    if subset
        DataWithinLimits = checkDataLimits( ...
            S(k).West, S(k).East, S(k).South, S(k).North, latlim, lonlim);
    else
        DataWithinLimits = 1;
    end
    
    if DataWithinLimits
        % Read the coordinate data.
        [S(k).Lat, S(k).Lon, S(k).BoundingBox] = readData(FileID, S(k));
    else
        % Move to the end of this data block.
        Offset = S(k).NumPoints*dataBlockLengthInBytes;
        fileseek(FileID, Offset, 'cof');
        
        % This element of the array S is not included in the output but
        % has been read. Set the counter back and reset the GSHHS_ID
        % number.
        S(k).GSHHS_ID = [];
        k = k-1;
    end
    
    % Update the FilePosition.
    FilePosition = ftell(FileID);  
end

% Remove empty elements.
index = cellfun(@isempty, {S.GSHHS_ID});
S(index) = [];

% Close the file.
fclose(FileID);

%--------------------------------------------------------------------------

function [readerFcn, extraAttributes] = getReaderFcn(FileID, version)
% Return a handle to a function that reads the header of the particular
% version of the GSHHS data.

v7ExtraAttributes = {'RiverLake', 'AreaFull', 'Container', 'Ancestor'};
switch version
    case 1
        readerFcn = @readV1Header;
        extraAttributes = {};
        
    case {2,3}
        readerFcn = @readV3Header;
        extraAttributes = {};
        
    case {4,5,6}
        readerFcn = @readV6Header;
        extraAttributes = {};
        
    case {7,8}
        readerFcn = @readV7Header;
        extraAttributes = v7ExtraAttributes;
        
    otherwise
        fname = fopen(FileID);
        latestQualifiedVersion = 8;
        warning('map:gshhs:unqualifiedVersion', ...
            ['The version number obtained from the file, ''%s'', is %d.', ...
            ' This version of the function %s is qualified to read data up through GSHHS version %d.', ...
            ' Attempting to continue to read the data from this version of the GSHHS dataset ...'], ...
            fname, version, upper(mfilename), latestQualifiedVersion)
        readerFcn = @readV7Header;
        extraAttributes = v7ExtraAttributes;
end

%--------------------------------------------------------------------------

function FileID = fileopen(filename)
% Open the filename and return a file identifier.

% Open the data file in binary mode, with big endian byte ordering.
FileID = fopen(filename,'rb','ieee-be');
if FileID == -1
    error('map:gshhs:invalidFileModes', ...
        'Unable to open file ''%s''.',filename)
end

%--------------------------------------------------------------------------

function FilePosition = fileseek(FileID, offset, origin)
% Wrap the fseek function to provide for current file position output
% and error checking.

status = fseek(FileID, offset, origin);
if status == -1
    error(ferror(FileID));
end
FilePosition = ftell(FileID);

%--------------------------------------------------------------------------

function [EOF, BOF] = getEofBofFilePositions(FileID)
% Get the end-of-file and beginning-of-file positions.

% Get the end of file position.
EOF = fileseek(FileID, 0, 'eof');
if EOF == 0
    fname = fopen(FileID);
    fclose(FileID);
    error('map:gshhs:noDataInFile', ...
        'No data in file:  ''%s''',  fname);
end

% Go back to the beginning of the file.
BOF = fileseek(FileID, 0, 'bof');

%--------------------------------------------------------------------------

function version = extractVersion(FileID)
% Extract the data version from the header and return it as a scalar
% double.

% Read the version number from the file.
version = readVersion(FileID);

% If the version number is not found in the file ([]), then assume
% that it is 1.
if isempty(version)
    version = 1;
end

%--------------------------------------------------------------------------

function numberOfHeaderBytes = getNumberOfHeaderBytes( ...
    FileID, readHeaderFcn, S, isWDB)
% Obtain number of bytes in the header.

[~, numberOfHeaderBytes] = readAttributes(FileID, readHeaderFcn, S, isWDB);
fileseek(FileID, -numberOfHeaderBytes, 'cof');

%--------------------------------------------------------------------------

function [S, numHeaderBytes] = readAttributes( ...
    FileID, readHeaderFcn, S, isWDB)
% Read the GSHHS attributes from the header records.

[header, numHeaderBytes] = readHeaderFcn(FileID);
validateHeader(header, FileID)
S = copyHeader(header, S, isWDB);

%--------------------------------------------------------------------------

function validateHeader(header, FileID)
% Validate the header.

if header.n <= 0 || header.n*2 > intmax('int32') % lat,lon pairs
    errorMsg = sprintf( ...
        'Invalid number of polygon points: %g.', header.n);
    attribError(FileID, 'invalidNumPoints', errorMsg)
end

if header.id < 0
    errorMsg = sprintf( ...
        'Invalid GSHHS ID: %g.', header.id);
    attribError(FileID, 'invalidId', errorMsg)
end

degreeCodeFactor = 1.0E06;   % Multiply by this to get back to degrees

limitSouth =  -91 * degreeCodeFactor;
limitNorth =   91 * degreeCodeFactor;
limitWest  = -181 * degreeCodeFactor;
limitEast  =  361 * degreeCodeFactor;

if header.south < limitSouth || header.south > limitNorth ...
        || header.north < limitSouth || header.north > limitNorth ...
        || header.east < limitWest || header.east > limitEast ...
        || header.west < limitWest || header.west > limitEast
    errorMsg = 'Invalid coordinate limits.';
    attribError(FileID, 'invalidCoordinateLimits', errorMsg)
end

if header.gcross > 1 || header.gcross < 0
    errorMsg = sprintf( ...
        'Invalid Greenwich crossing flag: %g.', header.gcross);
    attribError(FileID, 'invalidGreenwichCrossFlag', errorMsg)
end

if header.level < 0
    errorMsg = sprintf( ...
        'Invalid Level parameter: %g.', header.level);
    attribError(FileID, 'invalidLevelParameter', errorMsg)
end

if ~any(header.source == [0 1])
    errorMsg = sprintf( ...
        'Invalid Source parameter: %g.', header.source);
    attribError(FileID, 'invalidSourceParameter', errorMsg)
end

%--------------------------------------------------------------------------

function S = copyHeader(header, S, isWDB)
% Copy the header scalar structure to the scalar structure S.

if header.source == 0   % 0 = CIA WDBII, 1 = WVS
    source = 'WDBII';
else
    source = 'WVS';
end

switch header.level
    case 1
        levelString = 'land';
    case 2
        levelString = 'lake';
    case 3
        levelString = 'island_in_lake';
    case 4
        levelString = 'pond_in_island_in_lake';
    otherwise
        levelString = '';
end

if isWDB
    levelString = '';
    S.Geometry = 'Line';
end

areaDecodeFactor    = .1;      % Multiply by this to get back to km^2
degreeDecodeFactor = 1.0E-06;  % Multiply by this to get back to degrees

% Write the data to the output geographic data structure.
S.South = header.south * degreeDecodeFactor;
S.North = header.north * degreeDecodeFactor;
S.West =  header.west * degreeDecodeFactor;
S.East =  header.east * degreeDecodeFactor;
S.Area =  header.area * areaDecodeFactor;
S.Level = header.level;
S.LevelString     = levelString;
S.NumPoints       = header.n;
S.FormatVersion   = header.version;
S.Source          = source;
S.CrossGreenwich  = header.gcross;
S.GSHHS_ID        = header.id;

if isfield(S, 'RiverLake')
    S.RiverLake = header.river;
    S.AreaFull  = header.area_full;
    S.Container = header.container;
    S.Ancestor  = header.ancestor;
end
    
%--------------------------------------------------------------------------

function version = readVersion(FileID)
% Read the version number from the file opened by FileID. After obtaining
% the version number return the file pointer to the current position.

% Obtain current start position.
startFilePosition = ftell(FileID);

% Read the first 32 bytes from the header
H32 = read32ByteHeader(FileID);

% Check the second byte of the third 4-byte word to decide the format of
% data. In version 1.3 or earlier, this word will contain only the level
% value which ranges from 1 through 4. So the second byte will be all
% zeros. In versions 1.4 and later, this word contains four different
% values, each encoded in one 8-bit byte. The second byte contains the
% version number, which is always non-zero.
if H32(3) < 0
    attribError(FileID, 'inputsNotNonnegative', ...
    ['The header record is invalid. ', ...
    'Expecting a non-negative value for the version number.'])
end
headerVersion = bin2dec(dec2bin(bitget(H32(3),16:-1:9))');

if headerVersion == 0
    % For earlier GSHHS versions than 1.3 including 1.3
    % read the next int32 to obtain the version number.
    [version, count] = fread(FileID, 1, 'int32');
    if count ~= 1
        errorMsg = sprintf( ...
            ['Expecting to read 1 int32 header values.', ...
            ' Instead read %d int32 values.'], count);
        attribError(FileID, 'invalidInt32HeaderCount', errorMsg)
    end
    
    % Check the version number for earlier than 1.3
    if  version ~= 3
        version = [];
    else
        version = double(version);
    end
else
    % For GSHHS versions later than 1.3
    version = double(headerVersion);
end

% Return to startFilePosition.
endFilePosition = ftell(FileID);
offset = endFilePosition - startFilePosition;
fileseek(FileID, -offset, 'cof');

%--------------------------------------------------------------------------

function [header, numHeaderBytes] = readV1Header(FileID)

% Read the first 32 bytes from the header.
[H32, numHeaderBytes] = read32ByteHeader(FileID);

% Extract V1 header.
[header, numHeaderBytes] = extractV1Header(H32, FileID, numHeaderBytes);

%--------------------------------------------------------------------------

function [header, numHeaderBytes] = readV3Header(FileID)

% Read the first 32 bytes from the header.
[H32, numHeaderBytes] = read32ByteHeader(FileID);

% For versions older than V1 and including V1.3
[header, numHeaderBytes] = extractV3Header(H32, FileID, numHeaderBytes);

%--------------------------------------------------------------------------

function [header, numHeaderBytes] = readV6Header(FileID)

% Read the first 32 bytes from the header.
[H32, numHeaderBytes] = read32ByteHeader(FileID);

% For newer GSHHS versions than 1.3
[header, numHeaderBytes] = extractV4Header(H32, numHeaderBytes );

%--------------------------------------------------------------------------

function [header, numHeaderBytes] = readV7Header(FileID)

% Read the first 44 bytes from the header.
[H44, numHeaderBytes] = read44ByteHeader(FileID);

% For version 7 and higher.
[header, numHeaderBytes] = extractV7Header(H44, numHeaderBytes );

%--------------------------------------------------------------------------

function [H4, numHeaderBytes] = read4ByteHeader(FileID)
% Read 4 bytes from the header.

% Read INT16 header.
[H4, count] = fread(FileID,2,'int16');
if count ~= 2
    errorMsg = sprintf( ...
        ['Expecting to read 2 int16 header values.', ...
        ' Instead read %d int16 values.'], count);
    attribError(FileID, 'invalidInt16HeaderCount', errorMsg)
end
numHeaderBytes = count*2; % 2 bytes per int16

%--------------------------------------------------------------------------

function [H32,  numHeaderBytes] = read32ByteHeader(FileID)
% Read the first 32 bytes from the header.

[H32,count] = fread(FileID,8,'int32');
if count ~= 8
    errorMsg = sprintf( ...
        ['Expecting to read 8 int32 header values.', ...
        ' Instead read %d int32 values.'], count);
    attribError(FileID, 'invalidInt32HeaderCount', errorMsg)
end
numHeaderBytes = count*4; % 4 bytes per int32

%--------------------------------------------------------------------------

function [H44,  numHeaderBytes] = read44ByteHeader(FileID)
% Read the first 40 bytes from the header.

[H44,count] = fread(FileID,11,'int32');
if count ~= 11
    errorMsg = sprintf( ...
        ['Expecting to read 11 int32 header values.', ...
        ' Instead read %d int32 values.'], count);
    attribError(FileID, 'invalidInt44HeaderCount', errorMsg)
end
numHeaderBytes = count*4; % 4 bytes per int32

%--------------------------------------------------------------------------

function [header, numHeaderBytes] = extractV1Header(H32, FileID, numHeaderBytes)
% Extract version 1 header information.

% Read INT16 header attributes.
[H4, numBytes] = read4ByteHeader(FileID);
numHeaderBytes = numBytes + numHeaderBytes;

% Assign attributes to fields.
header.id     = H32(1);
header.n      = H32(2);
header.level  = H32(3);
header.gcross = H4(1);
header.source = H4(2);
header.west   = H32(4);
header.east   = H32(5);
header.south  = H32(6);
header.north  = H32(7);
header.area   = H32(8);

% Version number is not in included in the header, set to [].
header.version = [];

%--------------------------------------------------------------------------

function [header, numHeaderBytes] = extractV3Header(H32, FileID, numHeaderBytes)
% Extract version 3 header information.

% Read the next int32
[version, count] = fread(FileID, 1, 'int32');
if count ~= 1
    errorMsg = sprintf( ...
        ['Expecting to read 1 int32 header values.', ...
        ' Instead read %d int32 values.'], count);
    attribError(FileID, 'invalidInt32HeaderCount', errorMsg)
end

% Assign version number.
header.version = double(version);
numHeaderBytes = count*4 + numHeaderBytes;

% Read INT16 header attributes.
[H4, numBytes] = read4ByteHeader(FileID);
numHeaderBytes = numBytes + numHeaderBytes;

% Assign attributes to fields
header.id     = H32(1);
header.n      = H32(2);
header.level  = H32(3);
header.gcross = H4(1);
header.source = H4(2);
header.west   = H32(4);
header.east   = H32(5);
header.south  = H32(6);
header.north  = H32(7);
header.area   = H32(8);

%--------------------------------------------------------------------------

function [header, numHeaderBytes] = extractV4Header(H32, numHeaderBytes)
% Extract version 4 header information.

% Assign attributes to fields
header.id     = H32(1);
header.n      = H32(2);
flag = H32(3);
header.level   = bin2dec(dec2bin(bitget(flag, 8:-1:1))');
header.version = bin2dec(dec2bin(bitget(flag,16:-1:9))');
header.gcross  = bin2dec(dec2bin(bitget(flag,24:-1:17))');
header.source  = bin2dec(dec2bin(bitget(flag,32:-1:25))');
header.west   = H32(4);
header.east   = H32(5);
header.south  = H32(6);
header.north  = H32(7);
header.area   = H32(8);

%--------------------------------------------------------------------------

function [header, numHeaderBytes] = extractV7Header(H44, numHeaderBytes)
% Extract version 7 header information.

% Assign attributes to fields
header.id     = H44(1);
header.n      = H44(2);
flag = H44(3);
header.level   = bin2dec(dec2bin(bitget(flag, 8:-1:1))');
header.version = bin2dec(dec2bin(bitget(flag,16:-1:9))');
header.gcross  = bin2dec(dec2bin(bitget(flag,24:-1:17))');
header.source  = bin2dec(dec2bin(bitget(flag,24))');
header.river   = bin2dec(dec2bin(bitget(flag,32:-1:26))');

header.west   = H44(4);
header.east   = H44(5);
header.south  = H44(6);
header.north  = H44(7);
header.area   = H44(8);
header.area_full = H44(9);
header.container = H44(10);
header.ancestor  = H44(11);

%--------------------------------------------------------------------------

function [lat, lon, bbox] = readData(FileID, A)
% Read the latitude, longitude data and compute a bounding box.

[Data,count] = fread(FileID,[2,A.NumPoints],'int32');
if count ~= 2*A.NumPoints
    fclose(FileID);
    error('map:gshhs:dataReadError', ...
        ['Expecting to read %d latitude, longitude values. ', ...
        'Instead read only %d latitude, longitude values.'], ...
        2*A.NumPoints, count);
end

% Convert to degrees.
degreesPerMicroDegree = 1.0E-06;
lat = degreesPerMicroDegree*Data(2,:);
lon = degreesPerMicroDegree*Data(1,:);

% Fix longitude wrapping.
lon = radtodeg(unwrap(degtorad(lon)));
if any(lon > 195)
    lon = lon - 360;
end
if any(lon < -180)
    lon = lon + 360;
end

% Convert to clockwise ordering.
[lon, lat] = poly2cw(lon,lat);

% Fix Antarctica.
if A.South == -90
    if lon(1) == lon(end)
        lon(end) = [];
        lat(end) = [];
    end
    indexEast = find(lon < 180);
    indexWest = find(lon > 180);
    index = [fliplr(indexWest) fliplr(indexEast)];
    lon = lon(index);
    lat = lat(index);
    lon(lon > 180) = lon(lon > 180) - 360;
end

% Terminate with NaN.
lon(end+1) = NaN;
lat(end+1) = NaN;

% Compute the bounding box.
bbox = [ min(lon) min(lat); max(lon) max(lat)];

%--------------------------------------------------------------------------

function attribError(FileID, errorId, errorMsgPart2)
% Output common error message for attribute reading.

fname = fopen(FileID);
fclose(FileID);
error(['map:gshhs:' errorId], ...
    'Invalid attribute header in file "%s".\n%s', fname, errorMsgPart2)

%--------------------------------------------------------------------------

function S = initGeoStruct(numberOfRecords, extraNames)
% Create a geostruct to hold GSHHS data.

coordinateNames = {'Geometry', 'BoundingBox', 'Lat', 'Lon'};
attributeNames  = {'South', 'North', 'West', 'East', 'Area', ...
    'Level', 'LevelString', 'NumPoints', 'FormatVersion', ...
    'Source', 'CrossGreenwich', 'GSHHS_ID'};
names = [coordinateNames, attributeNames];

coordinateValues = {'Polygon', [], [], []};
attributeValues =  {[], [], [], [], [], '', '', [], '', '', [], []};
values = [coordinateValues, attributeValues];

if exist('extraNames', 'var')
    names = [names extraNames];
    values = [values cell(1, numel(extraNames))];
end

S(1:numberOfRecords,1) = cell2struct(values, names, 2);

%--------------------------------------------------------------------------

function tf = checkDataLimits(west,east,south,north,latlim,lonlim)
% Return false if the data in the current data block falls entirely outside
% the input limits (latlim & lonlim), and true if any part of the data
% falls inside the input limits.

full = abs(east - west) == 360;
if full
    longitudesOverlap = true;
else
    west = wrapTo180(west);
    east = wrapTo180(east);
    east(east < west) = east(east < west) + 360;
    longitudesOverlap = ~((lonlim(1) >=  east) | (lonlim(2) <=  west));
end

latitudesOverlap  = ~((latlim(1) >= north) | (latlim(2) <= south));
tf = latitudesOverlap & longitudesOverlap;

%--------------------------------------------------------------------------

function [extractindx, npolypts] = inlimitpolys(ifilename,ifilelength,latlim,lonlim)
% Check number of polygons in the file.

%total bytes in index file/( bits per number * number of numbers per record / bits/byte)
npoly = ifilelength/(32*5/8);

extractindx = [];
blocksize = 2000;
startblock = 0;
nrows  = npoly;
ncols = 5; % npts, latlim,lonlim

% Read number of points in each polygon from the index file.
readcols = [1 1];
readrows = [1 npoly];
npolypts = readmtx(ifilename,nrows,ncols,'int32',readrows,readcols,'ieee-be');

% Identify polygons in latlim. Do this in blocks to reduce memory
% requirements.
readcols = [2 5];
while 1
    
    readrows = [startblock*blocksize+1 min((startblock+1)*blocksize,npoly)];
    bbox = readmtx(ifilename,nrows,ncols,'int32',readrows,readcols,'ieee-be');
    bbox = bbox * 1.0E-06; % degrees (west east south north)
    
    % identify polygons that fall within the limits
    extractindx = ...
        [extractindx; ...
        (startblock*blocksize + ...
        find(checkDataLimits( bbox(:,1),bbox(:,2),bbox(:,3),bbox(:,4),...
        latlim,lonlim)) ) ...
        ];  %#ok<AGROW>
    
    if max(readrows) == npoly
        break
    end
    startblock = startblock+1;
end
