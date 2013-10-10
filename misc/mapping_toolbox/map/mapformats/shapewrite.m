function shapewrite(varargin)
%SHAPEWRITE Write geographic data structure to shapefile
%
%   SHAPEWRITE(S, FILENAME) writes a geographic data structure array S
%   (either a mapstruct or geostruct) to disk in shapefile format.  S
%   must be a valid geostruct or mapstruct, which implies the following
%   restrictions on its attribute fields:
%
%   * Each attribute field value must be either a real, finite, scalar
%     double or a character string.
%
%   * The type of a given attribute must be consistent across all features.
%
%   If S is a geostruct, then the values in its 'Lat' and 'Lon' fields
%   are written out as 'Y' and 'X' coordinates, respectively, in the
%   shapefile.
%
%   FILENAME must be a character string specifying the output
%   file name and location.  If an extension is included, it must be
%   '.shp' or '.SHP'. SHAPEWRITE creates three output files,
%
%                        [BASENAME '.shp']
%                        [BASENAME '.shx']
%                        [BASENAME '.dbf']
%
%   where BASENAME is FILENAME without its extension.
%
%   If a given attribute is integer-valued for all features, then it is
%   written to the [BASENAME '.dbf'] file as an integer.  If an attribute
%   is a non-integer for any feature, then it is written as a fixed point
%   decimal value with six digits to the right of the decimal place.
%
%   SHAPEWRITE(S, FILENAME, 'DbfSpec', DBFSPEC) writes a shapefile in which
%   the content and layout of the DBF file is controlled by a DBF
%   specification, indicated here by the parameter value DBFSPEC.  A DBF
%   specification is a scalar MATLAB structure with one field for each
%   feature attribute to be included in the output shapefile.  To include
%   an attribute in the output, make sure to provide a field in DBFSPEC
%   with a fieldname identical to the attribute name (the corresponding
%   fieldname in S), and assign to that field a scalar structure with the
%   following four fields:
% 
%     FieldName -- The field name to be used in the file
%
%     FieldType -- The field type to be used in the file ('N' or 'C')
%
%     FieldLength -- The field length in the file, in bytes
%
%     FieldDecimalCount -- For numeric fields, the number of digits to the
%                          right of the decimal place
%
%   When a DBF spec is provided, a given attribute will be included in the
%   output file only if it matches the name of a field in the spec.
%
%   The easiest way to construct a DBF spec is to call MAKEDBFSPEC, then
%   modify the output to remove attributes or change the FieldName,
%   FieldLength, or FieldDecimalCount for one or more attributes.  See the
%   help for MAKEDBFSPEC for more details and an example.  
%
%   Example
%   -------
%   % Derive a shapefile from concord_roads.shp in which roads of CLASS 5
%   % and greater are omitted.  Note the use of the 'Selector' option in
%   % shaperead, together with an anonymous function, to read only the main
%   % roads from the original shapefile.
%   shapeinfo('concord_roads')  % 609 features
%   S = shaperead('concord_roads', 'Selector', ...
%                 {@(roadclass) roadclass < 4, 'CLASS'});
%   shapewrite(S, 'main_concord_roads.shp')
%   shapeinfo('main_concord_roads')  % 107 features
%
%   See also MAKEDBFSPEC, SHAPEINFO, SHAPEREAD, UPDATEGEOSTRUCT.

% Copyright 2003-2009 The MathWorks, Inc.
% $Revision: 1.1.10.9 $  $Date: 2009/09/03 05:17:53 $

[S, basename, dbfspec] = parseInputs(varargin{:});
[shapeType, boundingBox, index] = writeSHP(S,basename);
writeSHX(shapeType, boundingBox, index, basename);
if ~isempty(dbfspec)
    dbfwrite(S, basename, dbfspec)
end

%----------------------------------------------------------------------------
function [shapeType, boundingBox, index] = writeSHP(S,basename)
% Write the main (SHP) file.

% Open the SHP file.
fid = fopen([basename '.shp'],'w','ieee-be');
if fid < 0
    eid = sprintf('%s:%s:failedToOpenSHP',getcomp,mfilename);
    error(eid,'Unable to open %s for writing.',[basename '.shp'])
end

% Get the shape type and a handle to function to write individual records.
shapeType = getShapeType(S(1).Geometry);
switch(S(1).Geometry)
    case 'Point'
        writefcn = @shpWritePoint;
    case 'MultiPoint'
        writefcn = @shpWriteMultiPoint;
    case {'Line','PolyLine','Polygon'}
        writefcn = @shpWritePoly;
end

% Obtain coordinates from the 'X' and 'Y' fields, if present.  Otherwise
% use the 'Lon' and 'Lat' fields.  Update file length and bounding box.
if isfield(S,'X') && isfield(S,'Y')
    xField = 'X';
    yField = 'Y';
elseif isfield(S,'Lon') && isfield(S,'Lat')
    xField = 'Lon';
    yField = 'Lat';
else
    eid = sprintf('%s:%s:missingCoordinateField',getcomp,mfilename);
    msg = 'Geographic data must have coordinate fields: either X and Y or Lat and Lon.';
    error(eid,msg)
end

% Write 100 bytes of zeros to reserve room for the file header.
fwrite(fid, uint8(zeros(1,100)), 'uint8');

% Write an SHP record for each element in S.
% Accumulate an index array and bounding box.
headerLengthInWords = 50;
fileLengthInWords = headerLengthInWords;
boundingBox = [Inf Inf -Inf -Inf];
index = zeros(2,numel(S));
for k = 1:numel(S)
    x = S(k).(xField);
    y = S(k).(yField);
    [fileLengthInWords, boundingBox, index(:,k)] = ...
        shpWriteRecord(fid, writefcn, shapeType, x, y, k, fileLengthInWords, boundingBox);
end

% Back up to the beginning of the file and write the header into the first 100 bytes.
fseek(fid, 0, 'bof');
shpWriteHeader(fid,shapeType,fileLengthInWords,boundingBox);

% Close the SHP file.
fclose(fid);

%----------------------------------------------------------------------------
function [fileLengthInWords, boundingBox, index] = shpWriteRecord(fid,...
      writefcn, shapeType, x, y, recordNumber, fileLengthInWords, boundingBox)
% Write an SHP record for each structure element s.  Return index entry.

recordHeaderLengthInWords = 4;
offsetInWords = ftell(fid) / 2;
[contentLengthInWords, recordBoundingBox] ...
    = writefcn(fid, recordNumber, shapeType, x, y);
fileLengthInWords = fileLengthInWords ...
    + recordHeaderLengthInWords + contentLengthInWords;
boundingBox(1:2) = min([boundingBox(1:2); recordBoundingBox(1:2)]);
boundingBox(3:4) = max([boundingBox(3:4); recordBoundingBox(3:4)]);
index = [offsetInWords; contentLengthInWords];

%----------------------------------------------------------------------------
function writeSHX(shapeType, boundingBox, index, basename)
% Write the SHX file.

% Open the SHX file.
fid = fopen([basename '.shx'],'w','ieee-be');
if fid < 0
    eid = sprintf('%s:%s:failedToOpenSHX',getcomp,mfilename);
    error(eid,'Unable to open %s for writing.',[basename '.shx'])
end

% Write the header.
headerLengthInWords = 50;
fileLengthInWords = headerLengthInWords + 4 * size(index,2);
shpWriteHeader(fid,shapeType,fileLengthInWords,boundingBox);

% Write the index records.
fwrite(fid,int32(index),'int32','ieee-be');

% Close the SHX file.
fclose(fid);

%----------------------------------------------------------------------------
function shapeType = getShapeType(geometry)

switch(geometry)
    case 'Point'
        shapeType = 1;
    case 'MultiPoint'
        shapeType = 8;
    case {'Line','PolyLine'}
        shapeType = 3;
    case {'Polygon'}
        shapeType = 5;
    otherwise
        shapeType = 0;
end

%----------------------------------------------------------------------------
function shpWriteHeader(fid, shapeType, fileLengthInWords, boundingBox)

fileCode = 9994;
version  = 1000;

bytes0thru27 = int32([fileCode 0 0 0 0 0 fileLengthInWords]);
fwrite(fid, bytes0thru27, 'int32', 'ieee-be');
% Check that count = 7

bytes28thru35 = int32([version shapeType]);
fwrite(fid, bytes28thru35, 'int32', 'ieee-le');
% Check that count = 2

bytes36thru99 = [boundingBox 0 0 0 0];
fwrite(fid, bytes36thru99, 'double', 'ieee-le');
% Check that count = 8

%----------------------------------------------------------------------------
function [contentLengthInWords, boundingBox] ...
    = shpWritePoint(fid,recordNumber,shapeType,X,Y)
% Write a Point record. Return the record content length measured in
% 16-bit words and the bounding box for the current record.

% Assign a degenerate bounding box for use in calculating the bounding box
% for the overall file.
boundingBox = [X Y X Y]; 

% The content length (the length of the record excluding its header) is
% fixed at 20 bytes (10 words) for Point shapes.
contentLengthInWords = 10;

% Write the record header
fwrite(fid, int32([recordNumber contentLengthInWords]), 'int32', 'ieee-be');

% Write the shape type.
fwrite(fid, int32(shapeType), 'int32', 'ieee-le');

% Write the point coordinates.
fwrite(fid, [X Y], 'double', 'ieee-le');

%----------------------------------------------------------------------------
function [contentLengthInWords, boundingBox] ...
    = shpWriteMultiPoint(fid,recordNumber,shapeType,X,Y)
% Write a MultiPoint record. Return the record content length measured in
% 16-bit words and the bounding box for the current record.

% Calculate the 2-by-numPoints points array and the bounding box.
points = [X(:)'; Y(:)'];
boundingBox = [min(X) min(Y) max(X) max(Y)]; 

% Calculate the content length (the length of the record excluding its header).
contentLengthInWords = (40 + 8 * numel(points)) / 2;

% Write the record header
fwrite(fid, int32([recordNumber contentLengthInWords]), 'int32', 'ieee-be');

% Write the shape type.
fwrite(fid, int32(shapeType), 'int32', 'ieee-le');

% Write the bounding box.
fwrite(fid, boundingBox, 'double', 'ieee-le');

% Write the number of points.
fwrite(fid, size(points,2), 'int32', 'ieee-le');

% Write the point coordinates.
fwrite(fid, points, 'double', 'ieee-le');

%----------------------------------------------------------------------------
function [contentLengthInWords, boundingBox] ...
    = shpWritePoly(fid,recordNumber,shapeType,X,Y)
% Write a PolyLine or Polygon record. Return the record content length measured in
% 16-bit words and the bounding box for the current record.

% Calculate the 1-by-numParts parts index, the 2-by-numPoints points array,
% and the bounding box.
[parts, points] = constructPartsIndexAndPointsArray(X,Y);
boundingBox = [min(X) min(Y) max(X) max(Y)]; 

% Calculate the content length (the length of the record excluding its header).
contentLengthInWords = (44 + 4 * numel(parts) + 8 * numel(points)) / 2;

% Write the record header.
fwrite(fid, int32([recordNumber contentLengthInWords]), 'int32', 'ieee-be');

% Write the shape type.
fwrite(fid, int32(shapeType), 'int32', 'ieee-le');

% Write the bounding box.
fwrite(fid, boundingBox, 'double', 'ieee-le');

% Write the part and point counts, and the parts index.
fwrite(fid, int32([numel(parts) size(points,2) parts]), 'int32','ieee-le');

% Write the point coordinates.
fwrite(fid, points, 'double', 'ieee-le');

%----------------------------------------------------------------------------
function [parts, points] = constructPartsIndexAndPointsArray(X,Y)
% Calculate the zero-based index array PARTS giving the offset to the where
% the start of each non-NaN sequence in vector X will be after all the NaNs
% were removed from X.  Organize the non-NaN elements of X and Y into a
% 2-by-numPoints array POINTS with X-coordinates in row 1 and Y-coordinates
% in row 2.

% Check that isequal(isnan(X),isnan(Y)). If not, then throw an error here
% and catch it above so we have a chance to close open files.  Also require
% that X and Y each have at least one non-NaN element.

% Determine where the NaNs are and add a (virtual) terminating NaN.
nanplaces = [find(isnan(X(:))); 1+numel(X)];

% Compute the parts index. Note that in both the initializion of parts and
% in subsequent updates, we calculate the offset to the next part before
% determining if that part exists.  This is fine, except that at the end we
% end up with parts(end) containing the offset of a non-existent part, so
% we must remove it.
parts = 0;
if any(nanplaces < numel(X))
    indexOfLastNaN = 0;
    for k = 1:numel(nanplaces)
        partLength = nanplaces(k) - (indexOfLastNaN + 1);
        if partLength > 0
            parts(end+1) = parts(end) + partLength; %#ok<AGROW>
        end
        % If partLength is zero, then the current NaN was preceded
        % immediately by another NaN, and should hence be ignored.
        indexOfLastNaN = nanplaces(k);
    end
    parts(end) = [];
end

% Calculate the points array:
%    [X(1) X(2) ...;
%     Y(1) Y(2) ...]
% after removing NaNs.
x = X(:)'; x(isnan(x)) = [];
y = Y(:)'; y(isnan(y)) = [];
points = [x; y];

%----------------------------------------------------------------------------
function dbfwrite(S, basename, dbfspec)
% Write the DBF file.

% Open the DBF file.
fid = fopen([basename '.dbf'],'w','ieee-le');
if fid < 0
    eid = sprintf('%s:%s:failedToOpenDBF',getcomp,mfilename);
    error(eid,'Unable to open %s for writing.',[basename '.dbf'])
end

% Determine the record length, construct format string, and for each field
% in the dbfspec, determine its position among the fields of the geographic
% data structure array, S.
[S, recordLength, fmt, geostructFieldIndex, nanmask] ...
    = setupRecordLayout(dbfspec, S);

% Write the DBF file: header, data records, and EOF marker
dbfWriteTableFileHeader(fid, numel(S), dbfspec, recordLength);
for k = 1:numel(S)
    dbfWriteTableRecord(fid, S(k), geostructFieldIndex, fmt, nanmask);
end
dbfMarkEOF(fid);

% Close the DBF file.
fclose(fid);

%----------------------------------------------------------------------------
function [S, recordLength, fmt, geostructFieldIndex, nanmask] = ...
    setupRecordLayout(dbfspec, S)

% While iterating to find the field length and building up a format
% specification string one attribute at a time, also construct a cell array,
% nulldata, which contains the value NaN for each numeric-valued
% attribute and the value '' for each string-valued attribute.
% This function modifies the string-valued attributes in the input
% geostruct array, S. For each character type attribute, convert the
% Unicode characters to their native representation, using the default
% character encoding stream.

geostructFields = fields(S);
recordLength = 1;  % Account for leading space
dbfFieldSpecs = struct2cell(dbfspec);
fmt = cell(size(dbfFieldSpecs));
nulldata = fmt;
attributeNames = fields(dbfspec);
geostructFieldIndex = zeros(size(attributeNames));
for k = 1:numel(dbfFieldSpecs)
    f = dbfFieldSpecs{k};
    switch(f.FieldType)
        case 'N'
            if f.FieldDecimalCount == 0
                fmt{k} = sprintf('%s%dd', '%', f.FieldLength);
            else
                fmt{k} = sprintf('%s%d.%df', '%', f.FieldLength, f.FieldDecimalCount);
            end
            nulldata{k} = NaN;
        otherwise % 'C'
            fmt{k} = sprintf('%s%ds', '%-', f.FieldLength);
            nulldata{k} = '';
            
            % Convert the Unicode characters to the default encoding scheme
            attributeName = attributeNames{k};
            for n=1:numel(S)
               S(n).(attributeName) = ...
                   char(unicode2native(S(n).(attributeName)));
            end
    end
    recordLength = recordLength + f.FieldLength;
    geostructFieldIndex(k) = ...
        strmatch(attributeNames{k}, geostructFields, 'exact');
end

% Concatenate the format strings for the individual fields into a single
% format string for the entire record, including a leading space.
fmt = [' ', fmt{:}];

% nanmask is a string, with length equal to recordLength.  It is used
% in subfunction dbfWriteTableRecord to ensure that each instance of a
% NaN-valued attribute is represented by blanks (space characters) in the
% xBase (.DBF) file.  nanmask consists of substrings equal to 'NaN'
% interspersed among a background of space characters (' ').  There is
% one such substring for each numerical attribute.  Each 'NaN' substring
% located within nanmask at exactly the location where 'NaN' would be
% written by a call to sprintf(fmt,c) if c were a cell array with a
% numerical value of NaN for the corresponding attribute field.  nanmask
% might look something like this:  '     NaN  NaN        NaN    ', for
% example, assuming three numeric-valued attributes.
nanmask = sprintf(fmt, nulldata{:});

%----------------------------------------------------------------------------
function dbfWriteTableFileHeader(fid, numberOfRecords, dbfspec, recordLength)

version = 3;

dateVector = datevec(date);
year  = dateVector(1);
month = dateVector(2);
day   = dateVector(3);

dbfFieldSpecs = struct2cell(dbfspec);
nFields = numel(dbfFieldSpecs);
headerLength = 32 + 32 * nFields + 1;

% Bytes 0-3: dBASE version and today's date
fwrite(fid, [version year-1900 month day], 'uint8');

% Bytes 4-7: Number of records in the table
fwrite(fid, numberOfRecords, 'uint32');

% Bytes 8-9: Number of bytes in the header
fwrite(fid, headerLength, 'uint16');

% Bytes 10-11: Number of bytes per record
fwrite(fid, recordLength, 'uint16');

% Bytes 12-31: Reserved (zero-fill)
fwrite(fid, zeros(1,20), 'uint8');

% Bytes 32-n: Table field descriptors (n = 32 + 32 * nFields)
for k = 1:nFields
    dbfWriteTableFieldDescriptor(fid, dbfFieldSpecs{k});
end

% Byte n+1: Field terminator
fwrite(fid, hex2dec('0D') ,'uint8');

%----------------------------------------------------------------------------
function dbfWriteTableFieldDescriptor(fid, dbfFieldSpec)

% Bytes 0-10: Field name in ASCII
%   (truncate or zero-fill to precisely 11 bytes)
fwrite(fid, truncateOrFill(dbfFieldSpec.FieldName,11), 'uchar');

% Byte 11: Field type in ASCII ('C' or 'N')
fwrite(fid, dbfFieldSpec.FieldType, 'uchar');

% Bytes 12-15: Field data address (zero-fill)
fwrite(fid, zeros(1,4), 'uint8');

% Byte 16: Field length
fwrite(fid, dbfFieldSpec.FieldLength, 'uint8');

% Byte 17: Field decimal count
%   (number of digits to the right of the decimal place)
fwrite(fid, dbfFieldSpec.FieldDecimalCount, 'uint8');

% Bytes 18-31: Reserved or zero
%   ("work area ID, SET FIELDS flag", ".MDX field flag")
fwrite(fid, zeros(1,14), 'uint8');

%----------------------------------------------------------------------------
function str = truncateOrFill(str,len)
% Ensure a string containing precisely LEN characters,
% padded with zeros if necessary.

str(len+1:end) = [];
str(end+1:len) = 0;

%----------------------------------------------------------------------------
function dbfWriteTableRecord(fid, s, geostructFieldIndex, fmt, nanmask)
% Write a single record, including the preceding SPACE character.
% s is a scalar geostruct or mapstruct, geostructFieldIndex is an array of
% indices indicating which fields of s need to be written, and
% fmt is a format string with one element per field in
% geostructFieldIndex.  See the comment in subfunction setupRecordLayout
% for a definition of nanmask.

% Note that if s contains any NaN-valued numerical attributes, then the
% call to sprintf will result in the substring 'NaN' being written into
% the record string.  Any such substrings must be replaced with blanks
% (space characters) before writing the record to the xBase (.DBF) file.
% This is accomplished via logical indexing with the nanmask string.

s = struct2cell(s);
record = sprintf(fmt, s{geostructFieldIndex});
record(record == nanmask) = ' ';
fwrite(fid, record, 'uchar');

%----------------------------------------------------------------------------
function dbfMarkEOF(fid)
% Mark end of dBASE file.

eofMarker = hex2dec('1A');
fwrite(fid, eofMarker, 'uint8');

%--------------------------------------------------------------------------
function [S, basename, dbfspec] = parseInputs(varargin)

validParameterNames = {'DbfSpec'};
numRequiredArgs = 2;
checknargin(numRequiredArgs,...
    numRequiredArgs + 2*numel(validParameterNames), nargin, mfilename);

% S, a geostruct or mapstruct array, is the first of two required inputs.
% It's expensive to validate, so we wait until we've checked
% everything else (except that we'll make sure it's not empty).
S = varargin{1};
if isempty(S)
    eid = sprintf('%s:%s:emptyGeostruct', getcomp, mfilename);
    error(eid, 'S must be a non-empty geographic data structure array.')
end

% FILENAME is the second of two required inputs.
filename = varargin{2};
basename = validateFilename(filename);

% Identify and validate the parameter name-value pairs beginning with
% the third argument.
for k = 3:2:nargin
    parName = checkstrs(varargin{k}, validParameterNames, ...
                        mfilename, sprintf('PARAM%d',(k-1)/2), k);                    
    switch parName
      case 'DbfSpec'
        checkExistence(k, nargin, mfilename, 'vector of record numbers', parName);
        dbfspec = varargin{k+1};
        
      otherwise
        eid = sprintf('%s:%s:internalProblem',getcomp,mfilename);
        error(eid,'Internal problem: unrecognized parameter name: %s',parName);
    end
end

% If the user provided a DBF spec, validate it.  Otherwise, derive one. 
if exist('dbfspec','var')
    % Validate dbfspec, validate the attribute fields of S, and ensure
    % for consistency between dbfspec and S.
    dbfspec = validateDbfSpec(dbfspec, S);
else
    % Note: makedbfspec validates the attribute fields of S.
    dbfspec = makedbfspec(S);
end

% Validate the geometry fields of S.
validateGeometry(S)

%--------------------------------------------------------------------------
function checkExistence(position, nargs, fcnName, ...
                        propertyDescription, propertyName)
% Error if missing the property value following a property name.

if (position + 1 > nargs)
    eid = sprintf('%s:%s:missingParameterValue',getcomp,fcnName);
    error(eid,...
          'Expected %s to follow parameter name ''%s''.',...
          propertyDescription, propertyName);
end

%--------------------------------------------------------------------------
function basename = validateFilename(filename)

checkinput(filename, {'char'}, {'vector'}, mfilename, 'FILENAME', 2);
[pathstr, name, ext] = fileparts(filename);
if ~any(strcmpi(ext,{'','.shp'}))
    eid = sprintf('%s:%s:invalidExtension',getcomp,mfilename);
    error(eid,'Function %s expected filename\n  %s\nto have the extension: ''shp''.',...
        mfilename, filename)
end
basename = fullfile(pathstr,name);

%----------------------------------------------------------------------------
function dbfspec = validateDbfSpec(dbfspec, S)

% If dbfspec is empty, make sure it's an empty struct.
if isempty(dbfspec)
    dbfspec = struct([]);
    return  % No need to check anything else, including the attribute fields of S.
end

% Make sure that dbfspec is a structure.
if ~isstruct(dbfspec)
    eid = sprintf('%s:%s:dbfspecNotAStructure', getcomp, mfilename);
    error(eid, 'dbfspec must be a structure.')
end

% Make sure that dbfspec is a scalar (or empty).
if numel(dbfspec) > 1
    eid = sprintf('%s:%s:nonScalarDbfSpec', getcomp, mfilename);
    error(eid, 'dbfspec must be a scalar (or empty) structure.')
end

% Validate S, then make sure that dbfspec and S are mutually consistent.
defaultspec = makedbfspec(S);  % Validates attribute values in S
attributeNamesInS = fields(defaultspec);
attributeNamesInDbfSpec = fields(dbfspec);

% Check 1:  Every attribute in dbfspec must be present in S.
missingAttributes = setdiff(attributeNamesInDbfSpec,attributeNamesInS);
if ~isempty(missingAttributes)
    eid = sprintf('%s:%s:missingAttributes', getcomp, mfilename);
    error(eid, ['dbfspec specifies attributes that are missing', ...
        ' from the geographic data structure S.'])
end

% Check 2:  Field types in dbfspec must be consistent with S.
%   While in this loop, use the default to fill in any fields missing from
%   dbfspec and ensure that the field length is at least 2 for
%   character-valued attributes and 3 for numerical attributes.
minCharFieldLength = 2;
minNumericalFieldLength = 3;  % Large enough to hold 'NaN'
for k = 1:numel(attributeNamesInDbfSpec)
    attributeName = attributeNamesInDbfSpec{k};
    fSpec = dbfspec.(attributeName);
    fDefault = defaultspec.(attributeName);
    
    if ~isfield(fSpec,'FieldName')
        dbfspec.(attributeName).FieldName = fDefault.FieldName;
    end

    if ~isfield(fSpec,'FieldType')
        dbfspec.(attributeName).FieldType = fDefault.FieldType;
    elseif fSpec.FieldType ~= fDefault.FieldType
        eid = sprintf('%s:%s:extaneousAttributes', getcomp, mfilename);
        error(eid, ['Type mismatch between dbfspec and geographic data', ...
            ' structure S for attribute: %s.'], attributeName)
    end
    
    if ~isfield(fSpec,'FieldLength')
        dbfspec.(attributeName).FieldLength = fDefault.FieldLength;
    elseif fSpec.FieldType == 'C'
        dbfspec.(attributeName).FieldLength ...
            = max(minCharFieldLength, fSpec.FieldLength);
    else
        dbfspec.(attributeName).FieldLength ...
            = max(minNumericalFieldLength, fSpec.FieldLength);
    end
    
    if ~isfield(fSpec,'FieldDecimalCount')
        dbfspec.(attributeName).FieldDecimalCount = fDefault.FieldDecimalCount;
    end
end

%----------------------------------------------------------------------------
function validateGeometry(S)

% Make sure there's a geometry field
if ~isfield(S, 'Geometry')
    eid = sprintf('%s:%s:noGeometry',getcomp,mfilename);
    error(eid,...
        'Geographic data structure S must have a ''Geometry'' field.');
end

% Check the coordinate field names
hasXY     = (isfield(S,'X') && isfield(S,'Y'));
hasLatLon = (isfield(S,'Lon') && isfield(S,'Lat'));
if ~hasXY && ~hasLatLon
    eid = sprintf('%s:%s:missingCoordinateFields',getcomp,mfilename);
    error(eid,'Geographic data structure S must have coordinate fields:\n%s',...
          'either ''X'' and ''Y'' or ''Lat'' and ''Lon''.')
end    
if (hasXY && hasLatLon)...
        || (hasXY     && (isfield(S,'Lat') || isfield(S,'Lon')))...
        || (hasLatLon && (isfield(S,'X')   || isfield(S,'Y')))
    eid = sprintf('%s:%s:redundantCoordinateFields',getcomp,mfilename);
    error(eid,'Geographic data structure S must have coordinate fields:\n%s',...
          'either ''X'' and ''Y'' or ''Lat'' and ''Lon'', but not both.')
end    

S1Geometry = S(1).Geometry;
% Make sure the geometry field of the first feature is a string
if ~ischar(S1Geometry)
    eid = sprintf('%s:%s:nonStrGeometry',getcomp,mfilename);
    error(eid,...
        ['The ''Geometry'' field of geographic data structure S', ...
        ' must contain a string.']);
end

% Make sure the geometry of the first feature has a valid value
validGeometries = {'Point', 'MultiPoint', 'Line', 'PolyLine', 'Polygon'};
if ~any(strcmp(S1Geometry,validGeometries))
    eid = sprintf('%s:%s:invalidGeometryString',getcomp,mfilename);
    error(eid,'%s\n%s',...
        ['The ''Geometry'' field of geographic data structure S', ...
        ' must contain one of the following strings:  ''Point'',', ...
        ' ''MultiPoint'', ''Line'', ''PolyLine'', or ''Polygon''.'])
end

% Check the coordinate values and geometry on a per-feature basis
if hasXY
    cfield1 = 'X';
    cfield2 = 'Y';
else
    cfield1 = 'Lon';
    cfield2 = 'Lat';
end
for k = 1:numel(S)
    if ~strcmp(S(k).Geometry,S(1).Geometry)
        eid = sprintf('%s:%s:inconsistentGeometry',getcomp,mfilename);
        error(eid,['All features in geographic data structure S', ...
            ' must have the same geometry.']);
    end
        
    c1 = S(k).(cfield1);
    c2 = S(k).(cfield2);
    
    % Make sure they're doubles
    if ~isa(c1,'double') || ~isa(c2,'double')
        eid = sprintf('%s:%s:nonDoubleCoordinate',getcomp,mfilename);
        error(eid,...
            'Coordinates in geographic data structure S must be doubles.')
    end
    
    % Make sure they match in length
    if ~isequal(size(c1),size(c2))
        eid = sprintf('%s:%s:coordinateSizeMismatch',getcomp,mfilename);
        error(eid,...
            ['Coordinate pairs in geographic data structure S', ...
            ' must match in size.'])
    end
    
    % Make sure they have matching NaNs, in any
    if isnan(c1) ~= isnan(c2)
        eid = sprintf('%s:%s:coordinateNanMismatch',getcomp,mfilename);
        error(eid,...
            ['Coordinate pairs in geographic data structure S', ...
            ' must have matching NaN locations.'])
    end
    
    % Make sure they're real and finite
    if ~all(isreal(c1)) || ~all(isreal(c2))...
            || ~all(isfinite(c1(~isnan(c1))))...
            || ~all(isfinite(c2(~isnan(c2))))
        eid = sprintf('%s:%s:coordinateNanMismatch',getcomp,mfilename);
        error(eid,...
            ['Coordinates in geographic data structure S', ...
            ' must have be real and finite.'])
    end
    
end

% For 'Point' geometry only, make sure all coordinate values are scalar
if strcmp(S(1).Geometry,'Point')
    for k = 1:numel(S)
        assert(isscalar(S(k).(cfield1)) && isscalar(S(k).(cfield2)), ...
            [getcomp ':' mfilename ':nonScalarPoint'], ...
            ['Coordinates in ''Point'' geographic data structure S', ...
            ' must be scalar.'])
    end
end
