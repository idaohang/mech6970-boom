function info = geotiffinfo(fileOrURL)
%GEOTIFFINFO Information about GeoTIFF file
%
%   INFO = GEOTIFFINFO(FILENAME) returns a structure whose fields contain
%   image properties and cartographic information about a GeoTIFF file.
%
%   FILENAME is a string that specifies the name of the GeoTIFF file.
%   FILENAME may include the directory name; otherwise, the file must be in
%   the current directory or in a directory on the MATLAB path.  If the
%   named file includes the extension '.TIF' or '.TIFF'  (either upper or
%   lower case), the extension may be omitted from  FILENAME.
%
%   If FILENAME is a file containing more than one GeoTIFF image, INFO is a
%   structure array with one element for each image in the file.  For
%   example, INFO(3) would contain information about the third image in the
%   file. If more than one image exists in the file it is assumed that each
%   image will have the same cartographic information and the same image
%   width and height.
%
%   INFO = GEOTIFFINFO(URL) reads the GeoTIFF image from an Internet URL.
%   The URL must include the protocol type (e.g., "http://").
%
%   The INFO structure contains the following fields:
%
%   Filename       A string containing the name of the file.
%
%   FileModDate    A string containing the modification date of the file.
%
%   FileSize       An integer indicating the size of the file in bytes.
%
%   Format         A string containing the file format, which should always
%                  be 'tiff'.
%
%   FormatVersion  A string or number specifying the file format version.
%
%   Height         An integer indicating the height of the image in pixels.
%
%   Width          An integer indicating the width of the image in pixels.
%
%   BitDepth       An integer indicating the number of bits per pixel.  
%
%   ColorType      A string indicating the type of image:
%                     'truecolor' for a truecolor (RGB) image
%                     'grayscale' for a grayscale image
%                     'indexed' for an indexed image
%
%   ModelType      A string indicating the type of coordinate system 
%                  used to georeference the image; either
%                     'ModelTypeProjected', 'ModelTypeGeographic' or ''.
%
%   PCS            A string describing the projected coordinate system.
%
%   Projection     A string describing the EPSG identifier for the 
%                  underlying projection method.
%
%   MapSys         A string indicating the map system; if applicable:
%                     'STATE_PLANE_27', 'STATE_PLANE_83', 
%                     'UTM_NORTH', 'UTM_SOUTH', or ''.
%
%   Zone           A double indicating the UTM or State Plane Zone number; 
%                  empty ([]) if not applicable or unknown.
%
%   CTProjection   A string containing the GeoTIFF identifier for the 
%                  underlying projection method.
%
%   ProjParm       A N-by-1 double containing projection parameter values.
%                  The identify of each element is specified by the 
%                  corresponding element of ProjParmId. Lengths are in
%                  meters, angles in decimal degrees.
%
%   ProjParmId     A N-by-1 cell array listing the projection parameter 
%                  identifier for each corresponding numerical element 
%                  ProjParm.  The possible values are:
%                     'ProjNatOriginLatGeoKey'   
%                     'ProjNatOriginLongGeoKey'
%                     'ProjFalseEastingGeoKey'   
%                     'ProjFalseNorthingGeoKey'
%                     'ProjFalseOriginLatGeoKey' 
%                     'ProjFalseOriginLongGeoKey'
%                     'ProjCenterLatGeoKey'      
%                     'ProjCenterLongGeoKey'
%                     'ProjAzimuthAngleGeoKey'   
%                     'ProjRectifiedGridAngleGeoKey'
%                     'ProjScaleAtNatOriginGeoKey'
%                     'ProjStdParallel1GeoKey'   
%                     'ProjStdParallel2GeoKey'
%
%   GCS            A string indicating the Geographic Coordinate System.
%
%   Datum          A string indicating the projection datum type, such as
%                  'North American Datum 1927' or 'North American Datum 1983'.
%
%   Ellipsoid      A string indicating the ellipsoid name as defined by
%                  the ellipsoid.csv EPSG file.
%
%   SemiMajor      A double indicating the length of the semi-major axis 
%                  of the ellipsoid, in meters.
%
%   SemiMinor      A double indicating the length of the semi-minor axis
%                  of the ellipsoid, in meters.
%
%   PM             A string indicating the prime meridian location, for
%                  example, 'Greenwich' or 'Paris'.
%
%   PmLongToGreenwich A double indicating the decimal degrees of longitude 
%                  between this prime meridian and Greenwich.  Prime 
%                  meridians to the west of Greenwich are negative.
%
%   UOMLength      A string indicating the units of length used in the
%                  projected coordinate system.
%
%   UOMLengthInMeters A double defining the UOMLength unit in meters.
%
%   UOMAngle       A string indicating the angular units used for 
%                  geographic coordinates.
%
%   UOMAngleInDegrees A double defining the UOMAngle unit in degrees.
%
%   TiePoints      A structure containing the image tiepoints.
%                  The structure contains these fields:
%
%                  ImagePoints  A structure containing the image row and 
%                               column coordinates of the tiepoints. 
%                               The structure contains these fields:
%                                  Row  A double array of size 1-by-N.
%                                  Col  A double array of size 1-by-N. 
%
%                  WorldPoints  A structure containing the x and y
%                               world coordinates of the tiepoints.
%                               The structure contains these fields:
%                                  X  A double array of size 1-by-N.
%                                  Y  A double array of size 1-by-N.
%                 
%   PixelScale     A 3-by-1 array that specifies the X,Y,Z pixel scale 
%                  values.
%
%   SpatialRef     A spatialref.MapRasterReference object if the model type 
%                  is 'ModelTypeProjected' or a
%                  spatialref.GeoRasterReference object if the model type
%                  is 'ModelTypeGeographic'. SpatialRef will be empty if
%                  the spatial referencing is ambiguously defined by the
%                  GeoTIFF file, or if the model type is
%                  'ModelTypeGeographic' and the geometric transformation
%                  type is 'affine'.
%
%   RefMatrix      The 3-by-2 referencing matrix. It must be unambiguously 
%                  defined by the GeoTIFF file, otherwise it will be empty.
%
%   BoundingBox    A 2-by-2 array that specifies the minimum (row 1) and 
%                  maximum (row 2) values for each dimension of the image 
%                  data in the GeoTIFF file.
%
%   CornerCoords   The CornerCoords structure contains the coordinates of 
%                  the outer corners of the GeoTIFF image.
%
%                  The CornerCoords structure contains six fields. Each is
%                  a 1-by-4 double array, or empty ([]), if unknown.  The
%                  arrays contain the coordinates of the outer corners of
%                  the corner pixels, starting from the (1,1) corner and 
%                  proceeding clockwise.
%
%                  X    Coordinates in the Projected Coordinate System.
%                       Equals Lon if the model type is 'ModelTypeGeographic'.
%
%                  Y    Coordinates in the Projected Coordinate System, 
%                       Equals Lat if the model type is 'ModelTypeGeographic'.
%
%                  Row  Row coordinates of the corner. 
%
%                  Col  Column coordinates of the corner.
%
%                  Lat  Latitudes of the corner. 
%
%                  Lon  Longitudes of the corner. 
%
%   GeoTIFFCodes   A structure containing raw numeric values for those 
%                  GeoTIFF fields which are encoded numerically in the
%                  file.  These raw values, converted to a string elsewhere
%                  in the INFO structure, are provided here for reference.
%
%                  The following fields are included:
%                     Model
%                     PCS
%                     GCS
%                     UOMLength
%                     UOMAngle
%                     Datum
%                     PM
%                     Ellipsoid
%                     ProjCode
%                     Projection
%                     CTProjection
%                     ProjParmId    
%                     MapSys
%
%                  Each is scalar except for ProjParmId which is a column
%                  vector.
%
%   GeoTIFFTags    A structure containing field names that match the 
%                  GeoTIFF tags found in the file. At least one GeoTIFF tag
%                  must be present in the file or an error is issued.
%                  The following fields may be included:
%
%                     ModelPixelScaleTag           1-by-3 double
%                     ModelTiepointTag             1-by-6 double
%                     ModelTransformationTag       1-by-16 double
%                     GeoKeyDirectoryTag           scalar structure
%                     GeoAsciiParamsTag            string
%                     GeoDoubleParamsTag           1-by-N double
%
%                  The GeoKeyDirectoryTag contains field names that match
%                  the names of the "GeoKeys". For more information about
%                  the "GeoKeys" refer to the GeoTIFF specification at
%                  http://www.remotesensing.org/geotiff/spec/geotiff6.html#6.2
%
%   ImageDescription A string describing the image; omitted if not included.
%
%   Example
%   -------
%   info = geotiffinfo('boston.tif')
%
%   See also GEOTIFFREAD, GEOTIFFWRITE, PROJFWD, PROJINV, PROJLIST.

% Copyright 1996-2010 The MathWorks, Inc.
% $Revision: 1.1.8.16.2.1 $  $Date: 2010/12/06 00:01:16 $

% Validate the input. If it's a filename, return the full pathname. If it's
% a URL, download the URL to a temporary file and set isURL to true.
[filename, isURL] = checkfilename( ...
    fileOrURL, {'tif', 'tiff'}, mfilename, 1, true);

% Create the information structure from contents in the file.
info = readinfo(filename);

% Delete temporary file from Internet download.
if (isURL)
    deleteDownload(filename);
end

%--------------------------------------------------------------------------

function info = readinfo(filename)
% Read the GeoTIFF info from the file and return the information structure.

% Obtain the TIFF information structure from the file.
tiff_info = tiffinfo(filename);

% Validate that at least one GeoTIFF tag exists in the file.
assert(~isempty(fieldnames(tiff_info(1).GeoTIFFTags)), ...
    'map:geotiffinfo:noGeoTiffTags', ...
    ['The file ''%s'' does not contain any GeoTIFF Tags. ',  ...
    'Use the function %s to obtain information about the file.'], ...
    filename, 'IMFINFO')
 
% Ensure that all images have the same Height and Width values.
assert(isconsistent(tiff_info.Height) && isconsistent(tiff_info.Width), ...
    'map:geotiffinfo:inconsistentImageSizes', ...
    'Multiple images exist in the file and their sizes are different.')

% Obtain the EPSG and PROJ directory.
[epsgDirName, projDirName] = feval(mapgate('getprojdirs'));

% Call the gtifinfo function to return the gtiff_info structure.
gtiff_info = gtifinfo(filename, epsgDirName, projDirName);

% Construct a new information structure and copy the input data to it.
info = constructInfo(tiff_info, gtiff_info);

%--------------------------------------------------------------------------

function tiff_info = tiffinfo(filename)
% Obtain subset of TIFF information from filename.

% Turn off TIFF library warnings for this function.
w = warning('off', 'MATLAB:tifflib:TIFFReadDirectory:libraryWarning');
wobj = onCleanup(@()warning(w));

% Open and validate the file.
% Close the TIFF file when the function terminates.
t = Tiff(filename);
tobj = onCleanup(@()close(t));
  
% Assign the tag names to retrieve from the file.
requiredTags = { ...
    'ImageLength', ...
    'ImageWidth', ...
    'BitsPerSample', ...
    'ColorMap', ...
    'Photometric'};
optionalTags = {'ImageDescription'};
geoTiffTags = getGeoTiffTagNames;
tiffTags = [requiredTags, optionalTags, geoTiffTags];

% Obtain the TIFF tag values from the Tiff object and return a structure
% with field names that match the names of the requested tags.
tags = getTags(t, tiffTags);

% Create the TIFF information structure.
tiff_info = tiffstruct(filename, tags,  optionalTags);

% Add the GeoTIFFTags field.
tiff_info = addGeoTIFFTagsField(tags, tiff_info, geoTiffTags);

%--------------------------------------------------------------------------

function geotiffNames = getGeoTiffTagNames
% Return a row cell vector of GeoTIFF tag names obtained from the Tiff class.

tiffNames = Tiff.getTagNames;
cindex1 = regexp(tiffNames, regexptranslate('wildcard', 'Model*Tag'));
cindex2 = regexp(tiffNames, regexptranslate('wildcard', 'Geo*Tag'));
index1 = ~cellfun(@isempty, cindex1); 
index2 = ~cellfun(@isempty, cindex2); 
geotiffNames = tiffNames(index1 | index2)';

%--------------------------------------------------------------------------

function tags = getTags(t, tiffTags)
% Obtain Tiff tag values from Tiff object, t.

% Obtain the number of TIFF directories.
numdirs = numberOfDirectories(t);

% Create a structure with the TIFF tags as field names.
tags = cell2struct(cell(size(tiffTags)), tiffTags, 2);
tags(numdirs) = tags(1);

% Obtain the TIFF tag values from the file.
for dirnum = 1:numdirs
    t.setDirectory(dirnum);
    for k=1:numel(tiffTags)
        tagname = tiffTags{k};
        try
            value = t.getTag(tagname);
            tags(dirnum).(tagname) = value;
        catch e
            id = 'MATLAB:tifflib';
            if ~strncmp(e.identifier, id, numel(id))
                rethrow(e);
            end
        end
    end
end

%--------------------------------------------------------------------------

function numdirs = numberOfDirectories(t)
% Obtain the number of valid TIFF directories.

numdirs = 1;
nofailures = true;
while ~t.lastDirectory && nofailures
    try
       t.nextDirectory;
       numdirs = t.currentDirectory;
    catch e
        if strcmp(e.identifier, 'MATLAB:tifflib:readDirectory:failure')
            % The next directory exists but cannot be read. numdirs will
            % not be updated and contains the number of valid directories
            % (one less than this directory number).
            nofailures = false;
        else
            % An unknown error occurred.
            rethrow(e)
        end
    end
end

% Return to the first directory.
t.setDirectory(1);

%--------------------------------------------------------------------------

function tiff_info = tiffstruct(filename, tags, optionalTags)
% Create the TIFF information structure. tags is a structure array with one
% element per TIFF directory. optionalTags is a cell array of strings.

% Create the information structure.
requiredFields = { ...
    'Filename', ...
    'FileModDate', ...
    'FileSize', ...
    'Format', ...
    'FormatVersion', ...
    'Height', ...
    'Width', ...
    'BitDepth', ...
    'ColorType'};

tiff_info = cell2struct(cell(size(requiredFields)), requiredFields, 2);
numdirs = numel(tags);
tiff_info(numdirs) = tiff_info;

% Obtain information about the file.
d = dir(filename);
fileModDate = d.date;
fileSize = d.bytes;
fileFormat = 'tif';

% Copy data to the tiff_info structure.
[tiff_info.Filename]      = deal(filename);
[tiff_info.FileModDate]   = deal(fileModDate);
[tiff_info.FileSize]      = deal(fileSize);
[tiff_info.Format]        = deal(fileFormat);
[tiff_info.FormatVersion] = deal([]);
[tiff_info.Height]        = deal(tags.ImageLength);
[tiff_info.Width]         = deal(tags.ImageWidth);

% Process BitDepth
for k=1:numel(tags)
    tiff_info(k).BitDepth = sum(tags(k).BitsPerSample);
end

% Process ColorType field.
for k = 1:numel(tags)
    if ~isempty(tags(k).ColorMap)
        tiff_info(k).ColorType = 'indexed';
    elseif tags(k).Photometric == Tiff.Photometric.MinIsWhite ...
            || tags(k).Photometric == Tiff.Photometric.MinIsBlack
        tiff_info(k).ColorType = 'grayscale';
    elseif tags(k).Photometric == Tiff.Photometric.RGB
        tiff_info(k).ColorType = 'truecolor';
    else
        tiff_info(k).ColorType = 'unknown';
    end
end

% Add the optionalTags if they are not empty.
for k=1:numel(optionalTags)
    tagname = optionalTags{k};
    if any(~cellfun(@isempty, {tags.(tagname)}))
        [tiff_info.(tagname)] = deal(tags.(tagname));
    end
end

%--------------------------------------------------------------------------

function tiff_info = addGeoTIFFTagsField(tags, tiff_info, geoTiffTags)
% Add the 'GeoTIFFTags' field to the info structure, tiff_info.

% Convert the GeoKeyDirectoryTag. Assume that all images have the same
% GeoKeyDirectoryTag information. (The MEX code, gtifinfo, parses the
% GeoTIFF tags on a file basis, not on an image basis).
if isfield(tags, 'GeoKeyDirectoryTag') ...
        && ~isempty(tags(1).GeoKeyDirectoryTag)
    tags(1).GeoKeyDirectoryTag = convertGeoKeys(tags(1));
end

% Add GeoTIFF tags to the tiff_info structure in the field 'GeoTIFFTags'.
% If there are no GeoTIFF tags in the file, the GeoTIFFTags structure is an
% empty structure. Assume that if the file contains multiple images, they
% have the same GeoTIFF tag information (consistent image size is enforced
% later). Translate the field name, 'GeoASCIIParamsTag' to
% 'GeoAsciiParamsTag' if present, based on the name in the GeoTIFF
% specification.
GeoTIFFTags = struct();
GeoTIFFTags(numel(tags)) = GeoTIFFTags;
for k=1:numel(geoTiffTags)
    tagname = geoTiffTags{k};
    if ~isempty(tags(1).(tagname))
        tagname_out = strrep(tagname, 'ASCII', 'Ascii');
        [GeoTIFFTags.(tagname_out)] = deal(tags(1).(tagname));
    end
end
[tiff_info.GeoTIFFTags] = deal(GeoTIFFTags);

%--------------------------------------------------------------------------

function GeoKeyDirectoryTag = convertGeoKeys(tags)
% Convert the GeoKeyDirectoryTag to a structure containing GeoKeys.

tag = tags.GeoKeyDirectoryTag;
assert(~isempty(tag), ...
    'map:geotiffinfo:GeoKeyDirectoryTagIsEmpty', ...
    'The %s cannot be empty.', 'GeoKeyDirectoryTag');

numKeys = tag(1,end);
assert(size(tag,2) == 4 && size(tag,1) >= numKeys+1, ...
    'map:geotiffinfo:invalidGeoKeyDirectoryTag', ...
    ['The %s does not have a valid size. Expected the size to be ', ...
    '%d-by-4 but the actual size is %d-by-%d.'], 'GeoKeyDirectoryTag', ...
    numKeys+1,  size(tag,1), size(tag,2));

% Construct the key ID to name map.
idMap = constructGeoKeyDirectoryMap;

% Obtain tag ID numbers for special data storage.
geoAscii  = Tiff.TagID.GeoASCIIParamsTag;
geoDouble = Tiff.TagID.GeoDoubleParamsTag;
asciiTag = 'GeoASCIIParamsTag';
if ~isfield(tags, asciiTag) || isempty(tags.(asciiTag))
    geoAsciiUnknown = 'GeoAsciiParamsTag is not in file.';
    haveAscii = false;
else
    haveAscii = true;
end
doubleTag = 'GeoDoubleParamsTag';
if ~isfield(tags, doubleTag) || isempty(tags.(doubleTag))
    geoDoubleUnknown = 'GeoDoubleParamsTag is not in file.';
    haveDouble = false;
else
    haveDouble = true;
end

% Translate the keys.
for k = 2:numKeys + 1
    id = tag(k,1);
    if idMap.isKey(id)
        key = idMap(id);
    else
        key = 'Unknown';
    end
    
    location = tag(k, 2);
    switch location
        case 0
            value = tag(k,end);
            
        case geoAscii
            if haveAscii
                value = tags.GeoASCIIParamsTag;
                count  = tag(k, 3) - 1; % Remove trailing |
                offset = tag(k, 4) + 1;
                startIndex = min(offset, numel(value));
                endIndex = min(offset+count-1, numel(value));               
                value = value(startIndex:endIndex);
            else
                value = geoAsciiUnknown;
            end
            
        case geoDouble
            if haveDouble
                value = tags.GeoDoubleParamsTag;
                count  = tag(k, 3);
                offset = tag(k, 4) + 1;
                startIndex = min(offset, numel(value));
                endIndex = min(offset+count-1, numel(value));    
                value = value(startIndex:endIndex);
            else
                value = geoDoubleUnknown;
            end
            
        otherwise
            value = 'Unknown';
    end
    GeoKeyDirectoryTag.(key) = value;
end

%--------------------------------------------------------------------------

function tf = isconsistent(varargin)
% True if there are two or more inputs that are numerically equal (see
% ISEQUAL) or if there are less than two inputs.

tf = (nargin < 2) || isequal(varargin{:});

%--------------------------------------------------------------------------

function info = constructInfo(tiff_info, gtiff_info)
% Create a new info structure and copy the input data to it.

% Set temporary height and width variables from the first element values.
height = tiff_info(1).Height;
width  = tiff_info(1).Width;

% Construct the spatial referencing object.
modelType = gtiff_info.ModelType;
tags = tiff_info(1).GeoTIFFTags(1);
SpatialRef = constructSpatialRef(tags, height, width, modelType);

% Construct the referencing matrix.
RefMatrix = constructRefMatrix(SpatialRef);

% Construct the TiePoints structure.
TiePoints = constructTiePoints(gtiff_info.TiePoints);

% Construct the BoundingBox structure.
BoundingBox = constructBoundingBox(SpatialRef);

% Construct the CornerCoords structure.
CornerCoords = constructCornerCoords(gtiff_info.CornerCoords, SpatialRef);

% Initialize a new structure.
info = initializeInfoStructure(tiff_info);

% Copy TIFF data.
info = copyTiffFields(info, tiff_info);

% Copy GeoTIFF data.
info = copyGeoTiffFields(info, gtiff_info);

% Copy other data.
[info.TiePoints]    = deal(TiePoints);
[info.RefMatrix]    = deal(RefMatrix);
[info.SpatialRef]   = deal(SpatialRef);
[info.BoundingBox]  = deal(BoundingBox);
[info.CornerCoords] = deal(CornerCoords);

%--------------------------------------------------------------------------

function SpatialRef = constructSpatialRef(Tags, height, width, modelType)
% Construct the SpatialRef object from TIFF tags.

useTiepoint = isfield(Tags, 'ModelTiepointTag') ...
    && isfield(Tags, 'ModelPixelScaleTag') ...
    && numel(Tags.ModelPixelScaleTag) >= 2 ...
    && numel(Tags.ModelTiepointTag) == 6;

useMatrix = isfield(Tags, 'ModelTransformationTag') ...
    && numel(Tags.ModelTransformationTag) >= 8;

fileHasSpatialData = (useTiepoint || useMatrix) ...
    && height ~= 0 && width ~= 0;

if fileHasSpatialData
    % The file contains valid spatial referencing data in the tags. Use the
    % values from the tags and the raster size to construct a spatial
    % referencing object.
    rasterSize = [height width];
    rasterInterpretation = getRasterInterpretation(Tags);
    SpatialRef = constructRasterReference( ...
        useMatrix, Tags, rasterSize, rasterInterpretation, modelType);
else
    % The file does not have valid referencing data in the tags or the
    % raster is invalid.
    SpatialRef = [];
end

%--------------------------------------------------------------------------

function SpatialRef = constructRasterReference( ...
    useMatrix, Tags, rasterSize, rasterInterpretation, modelType)
% Construct either a spatialref.MapRasterReference or
% spatialref.GeoRasterReference object based on the value of the string,
% modelType.

% Construct first corners and the Jacobian based on tag values.
if useMatrix
    transform = Tags.ModelTransformationTag;
    [x, y, J] = constructFromTransformMatrix(transform);
else
    tiepoint  = Tags.ModelTiepointTag;
    scale     = Tags.ModelPixelScaleTag;
    [x, y, J] = constructFromSingleTiepoint(tiepoint, scale);
end
       
switch modelType
    case 'ModelTypeGeographic'
        isRectilinear = (J(1,2) == 0 && J(2,1) == 0);
        if isRectilinear
            units = 'degrees';
            lat = y;
            lon = x;
            deltaLatNumerator = J(2,2);
            deltaLonNumerator = J(1,1);
            try
                SpatialRef = spatialref.GeoRasterReference( ...
                    rasterSize, rasterInterpretation, units, lat, lon, ...
                    deltaLatNumerator, 1, deltaLonNumerator, 1);
            catch e
                part1 = 'map:spatialref:';
                part2 = 'MATLAB:GeoRasterReference:';
                if strncmp(part1, e.identifier, numel(part1)) ...
                        || strncmp(part2, e.identifier, numel(part2))
                    SpatialRef = [];
                else
                    rethrow(e);
                end
            end
        else
            % A GeoRasterReference object cannot be created for a
            % (non-rectilinear) affine transformation.
            SpatialRef = [];
        end
        
    case 'ModelTypeProjected'
        units = '';
        try
            SpatialRef = spatialref.MapRasterReference(rasterSize, ...
                rasterInterpretation, units, x, y, J, [1 1; 1 1]);
        catch e
            part = 'map:spatialref:';
            if strncmp(part, e.identifier, numel(part))
                SpatialRef = [];
            else
                rethrow(e);
            end
        end
        
    otherwise
        SpatialRef = [];
end

%--------------------------------------------------------------------------

function rasterInterpretation = getRasterInterpretation(GeoTIFFTags)
% Obtain the raster interpretation from the GeoTIFF tags. From the tags, it
% will be 1 for 'RasterPixelIsArea' and 2 for 'RasterPixelIsPoint'. Convert
% the value to a string that is recognized by the RasterReference classes.

% If GTRasterTypeGeoKey ('RasterPixelIsArea' or 'RasterPixelIsPoint') is
% not specified, then default the value to 'RasterPixelIsArea'.
if ~isfield(GeoTIFFTags, 'GeoKeyDirectoryTag') ...
        || (isfield(GeoTIFFTags, 'GeoKeyDirectoryTag') ...
        && ~isfield(GeoTIFFTags.GeoKeyDirectoryTag, 'GTRasterTypeGeoKey'))
    rasterTypeGeoKey = 1; % RasterPixelIsArea
else
    rasterTypeGeoKey = GeoTIFFTags.GeoKeyDirectoryTag.GTRasterTypeGeoKey;
end

% Convert the value to a string.
if rasterTypeGeoKey == 1
    % 'RasterPixelIsArea'
    rasterInterpretation = 'cells';
else
    % 'RasterPixelIsPoint'
    rasterInterpretation = 'postings';
end

%--------------------------------------------------------------------------

function [x, y, J] = constructFromSingleTiepoint(tiepoint, scale)
% Construct the first corners and the Jacobian from a single tiepoint.

dx =  scale(1);
dy = -scale(2);
x = tiepoint(4) - dx * tiepoint(1);
y = tiepoint(5) - dy * tiepoint(2);
J = [dx 0; 0 dy];

%--------------------------------------------------------------------------

function [x, y, J] = constructFromTransformMatrix(transform)
% Construct the first corners and Jacobian from the transformation matrix.

x = transform(4);
y = transform(8);
J = transform([1 2; 5 6]);

%--------------------------------------------------------------------------
    
function CornerCoords = constructCornerCoords(Corners, SpatialRef)
% Construct the CornerCoords structure.

% Create a CornerCoords struct to hold the outer-edge corners.
CornerCoords = struct('X', [], ...
                      'Y', [], ...
                      'Row', [], ...
                      'Col', [], ...
                      'Lat', [], ...
                      'Lon', []);
                  
if ~isempty(SpatialRef)    
    % The Lat and Lon values will be calculated correctly from gtifinfo.
    cw = [1, 4, 3, 2];
    CornerCoords.Lat = Corners.Lat(cw)';
    CornerCoords.Lon = Corners.Lon(cw)';
    
    yi = SpatialRef.YLimIntrinsic([1 1 2 2]);
    xi = SpatialRef.XLimIntrinsic([1 2 2 1]);
    CornerCoords.Row = yi;
    CornerCoords.Col = xi;
    if strcmp('planar', SpatialRef.CoordinateSystemType)
        [xw, yw] = SpatialRef.intrinsicToWorld(xi, yi);
        CornerCoords.X = xw;
        CornerCoords.Y = yw;
    else
        CornerCoords.X = CornerCoords.Lon;
        CornerCoords.Y = CornerCoords.Lat;
    end   
end

%--------------------------------------------------------------------------

function BoundingBox = constructBoundingBox(SpatialRef)
% Construct the BoundingBox structure.
                  
if isempty(SpatialRef)   
   % The spatial referencing in the GeoTIFF file is not defined, or the 
   % image size is invalid, set BoundingBox to [].
   BoundingBox = [];   
else
   % Create the Bounding box from the spatial referencing object.
   if strcmp('planar', SpatialRef.CoordinateSystemType)
       x = SpatialRef.XLimWorld;
       y = SpatialRef.YLimWorld;
   else
       x = SpatialRef.Lonlim;
       y = SpatialRef.Latlim;
   end
   BoundingBox = [x' y'];
end

%--------------------------------------------------------------------------

function RefMatrix = constructRefMatrix(SpatialRef)
% Construct the referencing matrix, RefMatrix, from the spatial referencing
% object, SpatialRef.

if ~isempty(SpatialRef)
    W = SpatialRef.worldFileMatrix;
    
    C = [0  1  -1;...
         1  0  -1;...
         0  0   1];
    
    RefMatrix = (W * C)';
else
    RefMatrix = [];
end

%--------------------------------------------------------------------------

function TiePoints = constructTiePoints(tiepoints)
% Construct the TiePoints structure based on the values in the GeoTIFF 
% Tiepoint tag, tiepoints.

% Create an empty TiePoints struct.
imageCoordinates = struct('Row', [], 'Col', []);
mapCoordinates   = struct('X',   [], 'Y',   []);
TiePoints = struct('ImagePoints', imageCoordinates, ...
                   'WorldPoints', mapCoordinates);  

% n is the number of coordinates per tie point 
% (one triplet each of image and world coordinates)
n = 6;
numTiePoints = numel(tiepoints)/n;

if numTiePoints >= 1
    % The gtiff_info.TiePoints field is an array of 
    % GeoTIFF TiePoint coordinate values in the form 
    % (... I, J, K, X, Y, Z, ...); where
    % I, J, K correspond to the image coordinates, and
    % X, Y, Z correspond to the world coordinates.
    %
    % Referenced in the GeoTIFF specification at:
    % http://www.remotesensing.org/geotiff/spec/geotiff2.6.html#2.6
    %
    % Create the TiePoints struct with each field name
    % and set the corresponding GTIFF TiePoints values. 
    % Convert to 1-based image-coordinates.

    imageCoordinates.Row = (tiepoints(2:n:end) + .5)';
    imageCoordinates.Col = (tiepoints(1:n:end) + .5)';
    % Z (unused) is defined:
    % imageCoordinates.Z = [tiepoints(3:n:end)]';
    TiePoints.ImagePoints = imageCoordinates;

    mapCoordinates.X = (tiepoints(4:n:end))';
    mapCoordinates.Y = (tiepoints(5:n:end))';
    % Z (unused) is defined:
    % mapCoordinateStruct.Z = [gtiff_info.TiePoints(6:n:end)]';
    TiePoints.WorldPoints = mapCoordinates;
end

%--------------------------------------------------------------------------

function info = copyTiffFields(info, tiff_info)
% Copy TIFF info fields to info structure.
              
% The tiff_info structure may contain multiple elements.
tiffFields = fieldnames(tiff_info);
for k=1:numel(tiffFields)
   fieldName = tiffFields{k};
    [info.(fieldName)] = deal(tiff_info.(fieldName));
end

% If the ImageDescription field is set, copy it to the info structure.
if isfield(tiff_info,'ImageDescription')
   [info.ImageDescription] = deal(tiff_info.ImageDescription);
end

%--------------------------------------------------------------------------

function info = copyGeoTiffFields(info, gtiff_info)
% Copy GeoTIFF info fields to info structure.

% Set Zone value to [] if undefined.
if gtiff_info.Zone == 32767
   [gtiff_info.Zone] = deal([]);
end

% Remove unused fields and adjust field names.
gtiff_info = rmfield(gtiff_info, {'CornerCoords'});

% Copy the GeoTIFF info fields. The gtiff_info structure is scalar.
geoTiffFields = fieldnames(gtiff_info);
for k=1:numel(geoTiffFields)
    fieldName = geoTiffFields{k};
    [info.(fieldName)] = deal(gtiff_info.(fieldName));
end

%--------------------------------------------------------------------------

function info = initializeInfoStructure(tiff_info)
% Initialize info structure and set it's size to the size of tiff_info.

info = struct( ...
    'Filename', '', ...
    'FileModDate', '', ...
    'FileSize', [], ...
    'Format', '', ...
    'FormatVersion', [], ...
    'Height', [], ...
    'Width', [], ...
    'BitDepth', [], ...
    'ColorType', '', ...
    'ModelType', '', ...
    'PCS', '', ...
    'Projection', '', ...
    'MapSys', '', ...
    'Zone', [], ...
    'CTProjection', '', ...
    'ProjParm', [], ...
    'ProjParmId', [], ...
    'GCS', '', ...
    'Datum', '', ...
    'Ellipsoid', '', ...
    'SemiMajor', [], ...
    'SemiMinor', [], ...
    'PM', '', ...
    'PMLongToGreenwich', [], ...
    'UOMLength', '', ...
    'UOMLengthInMeters', [], ...
    'UOMAngle', '', ...
    'UOMAngleInDegrees', [], ...
    'TiePoints', [], ...
    'PixelScale', [], ...
    'SpatialRef', [], ...
    'RefMatrix',  [], ...
    'BoundingBox', [], ...
    'CornerCoords', [], ...
    'GeoTIFFCodes', [], ...
    'GeoTIFFTags', []);

if isfield(tiff_info,'ImageDescription')
   info.ImageDescription = '';
end
    
info(numel(tiff_info)) = info;
