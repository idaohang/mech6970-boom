function kmlwrite(varargin)
%KMLWRITE Write geographic data to KML file
%
%   KMLWRITE(FILENAME, LAT, LON) writes the latitude and longitude points,
%   LAT and LON, to disk in KML format. FILENAME must be a character string
%   specifying the output file name and location.  If an extension is
%   included, it must be '.kml'.  LAT and LON are numeric vectors, 
%   specified in degrees. LAT must be in the range [-90, 90]. There is no
%   range constraint on LON; all longitudes are automatically wrapped to
%   the range [-180, 180], to adhere to the KML specification.
%
%   KMLWRITE(FILENAME, S) writes a point or multipoint geographic data
%   structure to disk in KML format.  The Geometry field of S must be
%   either 'Point' or 'Multipoint'.  S must be a geostruct; that is, it
%   must include 'Lat' and 'Lon' coordinate fields.  (If S includes 'X'
%   and 'Y' fields an error is issued).  The attribute fields of S are
%   displayed as a table in the description tag of the placemark for
%   each element of S, in the same order as they appear in S.
%
%   KMLWRITE(FILENAME, ADDRESS) specifies the location of a KML Placemark
%   via an ADDRESS string or string cell array. Each string represents an
%   unstructured address with city, state, and/or postal code. If ADDRESS
%   is a cell array, each cell represents a unique point.
%
%   KMLWRITE(..., PARAM1, VAL1, PARAM2, VAL2, ...) specifies
%   parameter-value pairs that set additional KML feature properties.
%   Parameter names can be abbreviated and are case-insensitive. 
%
%   The parameter-value pairs are listed below:
%
%     Name
%             A string or cell array of strings which specifies a name
%             displayed in the viewer as the label for the object. If the
%             value is a string, the name is applied to all objects. If the
%             value is a cell array, it must have the same size as LAT and
%             LON, S, or ADDRESS.
%
%     Description
%             A string, cell array of strings, or attribute spec, which
%             specifies the contents to be displayed in the feature's
%             description tag(s). The description appears in the
%             description balloon when the user clicks on either the
%             feature name in the Google Earth Places panel or clicks the
%             placemark icon in the viewer window. If the value is a
%             string, the description is applied to all objects. If the
%             value is a cell array, it must have the same size as LAT and
%             LON, S, or ADDRESS.  Use a cell array to customize
%             descriptive tags for different placemarks.
%
%             Description elements can be either plain text or marked up
%             with HTML. When it is plain text, Google Earth applies basic
%             formatting, replacing newlines with <br> and giving anchor
%             tags to all valid URLs for the World Wide Web. The URL
%             strings are converted to hyperlinks. This means that you do
%             not need to surround a URL with <a href> tags in order to
%             create a simple link. Examples of HTML tags recognized by
%             Google Earth are provided on its Web site,
%             http://earth.google.com.
%
%             When an attribute spec is provided, the attribute fields of S
%             are displayed as a table in the description tag of the
%             placemark for each element of S.  The attribute spec is
%             ignored with LAT and LON input.  The attribute spec controls:
%
%                * Which attributes are included in the table
%                * The name for the attribute
%                * The order in which attributes appear
%                * The formatting of numeric-valued attributes
%
%             The easiest way to construct an attribute spec is to call
%             makeattribspec, then modify the output to remove attributes
%             or change the Format field for one or more attributes. 
%
%             Note that the 'Lat' and 'Lon' fields of S are not considered 
%             to be attributes. If included in an attribute spec, they are
%             ignored.
%
%     Icon
%             A string or cell array of strings which specifies a custom
%             icon filename. If the value is a string, the value is applied
%             to all objects. If the value is a cell array, it must have
%             the same size as LAT and LON, S, or ADDRESS.  If the icon
%             filename is not in the current directory, or in a directory
%             on the MATLAB path, specify a full or relative pathname. The
%             string may be an Internet URL. The URL must include the
%             protocol type (e.g., "http://"). 
%
%     IconScale
%             A positive numeric scalar or array which specifies a scaling
%             factor for the icon. If the value is a scalar, the value is
%             applied to all objects. If the value is an array, it must
%             have the same size as LAT and LON, S, or ADDRESS. 
%
%   View the KML file with a Google Earth browser
%   -----------------------------------------------
%   A KML file may be displayed in a Google Earth browser. Google Earth
%   must be installed on the system.
%
%   On Windows platforms, display the KML file with:
%      winopen(filename)
%
%   For Unix and MAC users, display the KML file with:
%      cmd = 'googleearth ';
%      fullfilename = fullfile(pwd, filename);   
%      system([cmd fullfilename])
%
%   View the KML file with a Google Maps browser
%   ----------------------------------------------
%   You can view KML files at the Google Maps Web site in addition to
%   Google Earth. The file must be located on a web server that is
%   accessible from the Internet.  A private intranet server will not
%   suffice, because Google's server must be able to access the URL that
%   you provide.  A template for using Google Maps is listed below:
%
%      GMAPS_URL = 'http://maps.google.com/maps?q=';
%      KML_URL = 'http://<your web server and path to your KML file>';
%      web([GMAPS_URL KML_URL])
%   
%   Example 1
%   ---------
%   % Write a single point to a KML file.
%   % Add a description containing HTML, a name and an icon.
%   lat =  42.299827;
%   lon = -71.350273;
%   description = sprintf('%s<br>%s</br><br>%s</br>', ...
%      '3 Apple Hill Drive', 'Natick, MA. 01760', ...
%      'http://www.mathworks.com');
%   name = 'The MathWorks, Inc.';
%   iconDir = fullfile(matlabroot,'toolbox','matlab','icons');
%   iconFilename = fullfile(iconDir, 'matlabicon.gif');
%   filename = 'MathWorks.kml';
%   kmlwrite(filename, lat, lon, ...
%      'Description', description, 'Name', name, 'Icon', iconFilename);
%
%   Example 2
%   ---------
%   % Write the locations of major European cities to a KML file, including 
%   % the names of the cities, and remove the default description table.
%   latlim = [ 30; 75];
%   lonlim = [-25; 45];
%   cities = shaperead('worldcities.shp','UseGeoCoords', true, ...
%      'BoundingBox', [lonlim, latlim]);
%   filename = 'European_Cities.kml';
%   kmlwrite(filename, cities, 'Name', {cities.Name}, 'Description',{});
%
%   Example 3
%   ---------
%   % Write the locations of several Australian cities to a KML file, 
%   % using addresses.
%   address = {'Perth, Australia', ...
%              'Melbourne, Australia', ...
%              'Sydney, Australia'};
%   filename = 'Australian_Cities.kml';
%   kmlwrite(filename, address, 'Name', address);
%
%   Example 4
%   ---------
%   % Write the locations of the Boston placenames to a KML file.
%   S = shaperead('boston_placenames');
%   proj = geotiffinfo('boston.tif');
%   surveyFeetPerMeter = unitsratio('sf','meter');
%   for k=1:numel(S)
%      x =  surveyFeetPerMeter * S(k).X;
%      y =  surveyFeetPerMeter * S(k).Y;
%      [S(k).Lat, S(k).Lon] = projinv(proj, x, y);
%   end
%   filename = 'Boston_Placenames.kml';
%   kmlwrite(filename, S, 'Name', {S.NAME});
%
%   See also GEOSHOW, MAKEATTRIBSPEC, SHAPEREAD, SHAPEWRITE.

% Copyright 2007-2009 The MathWorks, Inc.
% $Revision: 1.1.6.8 $  $Date: 2010/01/19 02:55:18 $

% Verify the number of varargin inputs.
error(nargchk(2,inf,nargin,'struct'));

% Verify the filename.
filename = verifyFilename(varargin{1});
varargin(1) = [];

% Parse the parameter/value pair input.
[numDataArgs, options, userSupplied, dataArgs] = parseOptions(varargin{:});

% Create a KML document object.
kmlDoc = KML_Document;  

% Set the name of the KML document.
[~, docName] = fileparts(filename);
setName(kmlDoc, docName);

% Add the data to the document.
addData(kmlDoc, numDataArgs, dataArgs, options, userSupplied);
 
% Create a KML file object.
kmlFile = KML_File(filename);  

% Write the KML document to the file.
write(kmlFile, kmlDoc); 

%--------------------------------------------------------------------------

function numDataArgs = verifyNumberOfDataArgs(varargin)
% Verify the number of required data arguments.
%
% The command line is composed of two sets of parameters:
%    the required arguments, followed by parameter-value pairs.
%
% For the case of KMLWRITE, the number of required arguments is either one
% or two. A single required argument may be of the following type: 
%   struct: geostruct input
%   char  : address data 
%   cell  : address data
%
% For the case of two required arguments (lat and lon input), the type of
% argument is numeric.
%
% verifyNumberOfDataArgs calculates the number of required data arguments
% and verifies the command line syntax. VARARGIN is expected to contain at
% least one element. The FILENAME argument from the command line is
% expected to be removed prior to calling.

% Assign a logical row vector that has one more element than VARARGIN,
% which contains true when the corresponding element of VARARGIN is a
% string, false otherwise, and which ends in true which will help in the
% case where VARARGIN contains no strings.
stringPattern = [cellfun(@ischar, varargin) true];

% Determine the number of data arguments:
%
% Define a working row vector that looks like this: [1 1 2 3 ...
% numel(varargin)], and call it numDataArgKey.  Then index into it with the
% stringPattern row vector (both have length 1 + numel(varargin)) and keep
% the first element that results (because there may be more than one
% string) to define numDataArgs.  The following three cases cover all the
% possibilities:
%
% (1) varargin{1} is a string, which means it must contain an address,
% which is the one and only data argument.  The first element of
% stringPattern is true, and the first element returned by the logical
% indexing step is the first element of numDataArgKey, which is 1.
%
% (2) varargin contains at least one string, but it is not varargin{1}.  In
% this case, the number of data arguments is one less than the position of
% the first string in varargin. Due to the offset in the definition of
% numDataArgKey (the 2nd value is 1, the 3rd value is 2, etc.) the first
% element of the array returned by the indexing step will be the required
% value.
%
% (3) varargin contains no strings, so only the last element in
% stringPattern is true, which indexes the last element of numDataArgKey,
% which is numel(varargin). All the arguments are data arguments in this
% case.
numDataArgKey = [1, 1:numel(varargin)];
numDataArgs = numDataArgKey(stringPattern);
numDataArgs = numDataArgs(1);

% Set logicals for the conditions of one or two data arguments. For one
% data argument, the char case is previously defined by numDataArgs == 1.
% If there is not one data parameter, then there must be two.
haveOneDataArg = ...
    numDataArgs == 1 || iscell(varargin{1}) || isstruct(varargin{1});
haveTwoDataArgs = ~haveOneDataArg;

% The paramPos is the expected position of the first parameter-value pair.
% The function expects only one or two data arguments, so the paramPos is
% one plus the expected values.
paramPos = [2,3];
paramPos = paramPos([haveOneDataArg haveTwoDataArgs]);

% If the paramPos variable is not a string, then error.
if ~stringPattern(paramPos)

    % argPos is the position in the command line arguments for the paramPos
    % variable. Since FILENAME is removed from VARARGIN, it is one plus
    % paramPos. 
    argPos = num2ordinal(paramPos+1);
    error('map:kmlwrite:invalidType', ...
        ['Function KMLWRITE expected its %s input, PARAM,\n', ...
        'to be one of these types:\n\nchar\n\n', ...
        'Instead its type was %s'], argPos, class(varargin{paramPos}));
end
    
%--------------------------------------------------------------------------

function [numDataArgs, options, userSupplied, dataArgs] = ...
    parseOptions(varargin)
% Parse varargin for parameter/value pairs. Return the parameters as field
% names in options. userSupplied is a struct with parameters as field
% names. The value is set to true if the parameter is supplied in varargin.
% Remove the pairs and return the result in dataArgs. numDataArgs is the
% number of arguments to the first parameter-value pair.

% Verify the number of data arguments.
numDataArgs = verifyNumberOfDataArgs(varargin{:});

% Define dataArgs.
dataArgs = varargin(1:numDataArgs);

% Remove the data arguments from varargin to obtain the pvpairs.
varargin(1:numDataArgs) = [];

% Set the number of required elements for the options structure. If a
% single character address is specified, then the number of elements is 1;
% otherwise, the number of elements is equal to the numel of the first data
% argument.
haveCharAddressData = ischar(dataArgs{1});
if haveCharAddressData
    numElements = 1;
else
    numElements = numel(dataArgs{1});
end

% parameterNames is a cell array of valid parameter names.
parameterNames = {'Description','Name','Icon','IconScale'};

% validateFcns is a cell array of a validation functions for each
% parameterName.
validateFcns = { ...
    @(x)validateCellWrapper(x, @validateDescription, ...
       parameterNames{1}, numElements), ...
    @(x)validateCellWrapper(x, @validateStringCell, ...
       parameterNames{2}, numElements), ...
    @(x)validateCellWrapper(x, @validateFilenameCell, ...
       parameterNames{3}, numElements), ...
    @(x)validateNumericWrapper(x, @validatePositiveNumericArray, ...
       parameterNames{4}, numElements)};
 
% Parse the parameter/value pair input.
[opt, userSupplied, unmatched] = ...
    internal.map.parsepv(parameterNames, validateFcns, varargin{:});

% Check if varargin contained unmatched parameters.   
if ~isempty(unmatched)
    permissibleParams = ...
        [sprintf('''%s'', ', parameterNames{1:end-1}), ...
        'or ' '''' parameterNames{end} ''''];
    error('map:parsepv:pvpairCountMismatch', ...
        ['Function KMLWRITE expected the parameter, ''%s'', to be one of the following:\n' ...,
        permissibleParams '.'], unmatched{1});
end

% opt is a scalar struct. Expand to a struct array, options, matching in
% size to the number of elements of S or the number of elements of LAT and
% LON.
options(max(numElements,numel(opt)),1) = opt(1);

% Loop through each parameterName to convert options to a structure array.
for k=1:numel(parameterNames)
    % Set the parameterName from the parameterNames list.
    parameterName = parameterNames{k};
   
    % Convert the cell array of each parameter (opt.(parameterName)) to a
    % struct array, with one element in the field parameterName. Set the
    % field of the options struct to the output structure array.
    S = cellToStructArray(opt.(parameterName), parameterName, numElements);
    [options.(parameterName)] = S.(parameterName);
end

% If an attribSpec is supplied in the options structure, then set the Spec
% description field to it, otherwise set it to []. If an attribSpec is
% provided, then the userSupplied.Description field needs to be set to
% false, since the Description will be created using the supplied
% attribSpec.
if userSupplied.Description && isstruct(options(1).Description)
    options(1).Spec = options(1).Description;   
    options(1).Description = {' '};
    userSupplied.Description = false;
else
    options(1).Spec = [];
end

%--------------------------------------------------------------------------

function c = validateCellWrapper(c, validateFcn, parameter, numElements)
% Validate cell wrapper function provides a common interface to validate
% inputs that are required to be cell array.

% c needs to be a cell array.
if ~iscell(c)
    c = {c};
end

% c needs to be a row vector.
c = {c{:}}; %#ok<CCAT1>

% Validate the number of elements in c.
validateNumberOfCellElements(c, parameter, numElements);

% Execute the validation function.
c = validateFcn(c, parameter);

% Map any empty characters to space characters.
c = mapEmptyToSpaceChar(c);

%--------------------------------------------------------------------------

function c = validateNumericWrapper(c, validateFcn, parameter, numElements)
% Validate numeric wrapper function provides a common interface to validate
% inputs that are required to be numeric array.

if iscell(c)
    c = cell2mat(c);
end

% Execute the validation function.
c = validateFcn(c, parameter);

% All parameter values must be converted to a cell array.
c = num2cell(c(:)');

% Validate the number of elements in c.
validateNumberOfCellElements(c, parameter, numElements);

% Map any empty characters to space characters.
c = mapEmptyToSpaceChar(c);

%--------------------------------------------------------------------------

function c = validateDescription(c, parameter)
% Validate description input. c is a cell array that is validated to
% contain strings or struct input.

cIsCellArrayOfStringsOrStruct = ...
    all(cellfun(@(x)(ischar(x) || isstruct(x)), c));
assert(cIsCellArrayOfStringsOrStruct, ...
    'map:kmlwrite:expectedStringOrStructCell', ...
    ['Function KMLWRITE expected ''%s'' value to be class type\n', ...
    'string or struct.'], parameter);

%--------------------------------------------------------------------------

function  c = validateStringCell(c, parameter)
% Validate c to be a cell array of strings.

cIsCellArrayOfStrings = all(cellfun(@ischar, c));
assert(cIsCellArrayOfStrings, ...
    'map:kmlwrite:expectedStringCell', ...
    'Function KMLWRITE expected ''%s'' value to be class type string.', ...
    parameter);

%--------------------------------------------------------------------------

function  c = validatePositiveNumericArray(c, parameter)
% Validate c to be an array containing positive numeric values. 

cIsPositiveNumericArray = ...
    isnumeric(c) && all(~isinf(c(:))) && all(c(:) > 0);
assert(cIsPositiveNumericArray, ...
    'map:kmlwrite:expectedNumericArray', ...
    ['Function KMLWRITE expected the ''%s'' value\n', ...
    'to be a positive numeric class type.'], ...
    parameter);

%--------------------------------------------------------------------------

function  c = validateFilenameCell(c, parameter)
% Validate c to be a cell array containing filenames. The filenames are
% validated to be strings and to exist. A filename may contain a URL
% string containing ftp:// http:// or file://.

% Validate the input as all string.
c = validateStringCell(c, parameter);

% Verify that all files exist. Some files may be a URL string, in which
% case filesExist is set to false.
filesExist = logical(cellfun(@(x)exist(x,'file'), c));

% urlEntries is a logical array that is true for entries that contain a
% URL string.
urlEntries = isURL(c);

% filesAreValid is a logical array set to true for all entries that are
% valid.
filesAreValid = urlEntries | filesExist;

% invalidEntries is a cell array of entries that are invalid.
invalidEntries = c(~filesAreValid);

% invalidEntriesString is a char list of entries that are invalid.
invalidEntriesString = strtrim(sprintf(' %s', invalidEntries{:}));
invalidEntriesString = strrep(invalidEntriesString,' ', ', ');

assert(all(filesAreValid), ...
    'map:kmlwrite:invalidFilename', ...
    'Function KMLWRITE was unable to find icon file or files ''%s''.', ...
    invalidEntriesString);

% The files may be partial pathnames. Set all filenames to absolute path.
c = getAbsolutePath(c, filesExist);

%--------------------------------------------------------------------------

function filenames = getAbsolutePath(filenames, filesExist)
% Return the absolute path of each element in filenames. filesExist is a
% cell array that is true for each file that exists.

for k=1:numel(filenames)
    if filesExist(k)
        
        try 
           fid = fopen(filenames{k},'r');
           fullfilename = fopen(fid);
           fclose(fid);
        catch e
            error('map:kmlwrite:unableToOpen', ...
                'Function KMLWRITE was unable to open icon file ''%s''.', ...
                filenames{k});
        end
                
        if exist(fullfile(pwd,fullfilename),'file')
           fullfilename = fullfile(pwd, fullfilename);
        end
        filenames{k} = fullfilename;
    end
end
    
%--------------------------------------------------------------------------

function tf = isURL(filenames)
% Determine if a cell array of filenames contain a URL string. Return a
% logical array that is true for each element in filenames that contains a
% URL string.

urlFiles = strfind(filenames, '://');
tf = cellfun(@isempty, urlFiles);
tf = ~tf;

%--------------------------------------------------------------------------

function validateNumberOfCellElements(c, parameter, maxNumElements)
% Validate the number of elements in the c cell array.

validNumberOfCellElements ...
    = ismember(numel(c),[0, 1, maxNumElements]);
assert(validNumberOfCellElements, ...
    'map:kmlwrite:mismatchNumelOptions', ...
    ['Function KMLWRITE expected the number of value elements of ''%s''\n', ...
    'to be 1 or to match the number of elements, ''%d'', in the data input.'], ...
    parameter, maxNumElements);

%--------------------------------------------------------------------------

function c = mapEmptyToSpaceChar(c)
% Map empty values to a single space character. c is a cell array. empty
% values in the cell array are changed to ' '.
%
% XMLWRITE will not output empty tags correctly. For example, a value of ''
% for a tag name of 'description' will output as:
% <description/> 
% rather than:
% <description></description>

spaceChar = ' ';
c(isempty(c)) = {spaceChar};
emptyIndex = cellfun(@isempty, c);
c(emptyIndex) = {spaceChar};

%--------------------------------------------------------------------------

function S = cellToStructArray(c, name, numElements)
% Convert a cell array to a struct array. 
%
% The cell array, c, has either one or numElements elements, which must be
% validated prior to calling this function. The struct array, S, has
% numElements elements. S has a single field name, name.

% Expand c in size to numElements.
% To minimize space there is no need to expand if the array contains a spec
% structure.
if numel(c) == 1 && ~isstruct(c{1})
    c(1:numElements) = c(1);
elseif isstruct(c{1}) && numElements > 1
    c(2:numElements) = {' '};
end

% Convert c to a struct array with fieldname 'name'.
S = cell2struct(c, name, 1);

%--------------------------------------------------------------------------

function [lat, lon] = verifyCoordinates(lat, lon)
% Verify and validate the latitude, longitude coordinates.

% Validate numeric data.
isValidClass = @(x) isreal(x) && isnumeric(x);
assert(isValidClass(lat) && isValidClass(lon), ...
    'map:kmlwrite:nonNumericCoordinate', ...
    ['Function KMLWRITE expected the LAT and LON coordinate arrays\n', ...
    'to be real and numeric.']);

% The shape of the arrays do not matter.
lat = lat(:);
lon = lon(:);
    
% Validate finite.
assert(all(~isinf(lat)) && all(~isinf(lon)), ...
    'map:kmlwrite:isInfinite',  ...
    ['Function KMLWRITE expected the LAT and LON coordinate arrays\n', ...
    'to be finite.']);

% Validate size of arrays equal.
assert(numel(lat) == numel(lon), ...
    'map:kmlwrite:latLonSizeMismatch', ...
    ['Function KMLWRITE expected the LAT and LON coordinate arrays\n', ...
    'to contain the same number of elements.']);

% Validate arrays are not empty.
assert(~isempty(lat), ...
    'map:kmlwrite:latLonIsEmpty',  ...
    ['Function KMLWRITE expected the LAT and LON coordinate arrays\n', ...
    'to be nonempty.']);

% Validate NaN locations.
latNaNIndices = isnan(lat);
lonNaNIndices = isnan(lon);
assert(isequal(latNaNIndices, lonNaNIndices), ...
    'map:kmlwrite:invalidNaNIndices',  ...
    ['Function KMLWRITE expected the NaN locations of the\n', ...
    'LAT and LON coordinate arrays to match in size and location.']);

% Validate at least one non-NaN coordinate.
assert(~all(latNaNIndices), ...
    'map:kmlwrite:latLonIsNaN', ...
    ['Function KMLWRITE expected the LAT and LON coordinate arrays\n', ...
    'to contain at least one non-NaN value.']);

% Validate latitudes are in range.
allLatValuesInRange = ...
    all(-90 <= lat(~latNaNIndices)) &&  all(lat(~latNaNIndices) <= 90);
assert(allLatValuesInRange, ...
    'map:kmlwrite:latOutOfRange', ...
    ['Function KMLWRITE expected the LAT coordinate array\n', ...
    'to be in the range [-90, 90].']);

%--------------------------------------------------------------------------

function filename = verifyFilename(filename)
% Verify and validate the filename input.

filenamePos = 1;
checkinput(filename, {'char'}, {'vector','nonempty'}, mfilename, ...
    'FILENAME', filenamePos);

% Add .kml extension to filename if necessary.
kmlExt = '.kml';
[pathname, basename, ext] = fileparts(filename);
assert(isequal(lower(ext), kmlExt) || isempty(ext), ...
    'map:kmlwrite:invalidExtension', ...
    ['Function KMLWRITE expected the file extension of ''%s''\n', ...
    'to be ''%s'' or empty.'], filename, kmlExt);
filename = fullfile(pathname,[basename,kmlExt]);

%--------------------------------------------------------------------------

function addData(...
    kmlDoc, numDataArgs, dataArgs, options, userSupplied)
%Add data to the KML document based on number of data arguments.

% We know that numDataArgs can only be either 1, or 2. 
switch numDataArgs
    case 1
        validInput = ...
            isstruct(dataArgs{1}) || iscell(dataArgs{1}) || ischar(dataArgs{1});
        assert(validInput, ...
            'map:kmlwrite:invalidSecondInput', ...
            ['Function KMLWRITE expected its second input, S or ADDRESS,\n', ...
            'to be one of these types:\n\nstruct, char, or string cell array\n\n', ...
            'Instead its type was %s'], class(dataArgs{1}));

        if isstruct(dataArgs{1}) 
            addStruct(kmlDoc, dataArgs{:}, options, userSupplied);
        else
            addAddress(kmlDoc, dataArgs{:}, options);
        end
            
    case 2
        addPoints(kmlDoc, dataArgs{:}, options);
end

%--------------------------------------------------------------------------

function addAddress(kmlDoc, address, options)
% Add addresses to the KML document.

% address must be a cell array.
if ~iscell(address)
    address = {address};
end

% Verify that address is a cell array of strings.
isCellArrayOfStrings = all(cellfun(@ischar, address));
assert(isCellArrayOfStrings, ...
    'map:kmlwrite:expectedStringCell', ...
    ['Function KMLWRITE expected its second input, ADDRESS, \n', ...
    'to be a string or string cell array of strings.']);

% Append the addresses to the KML document.
appendAddress(kmlDoc, address, options);

%--------------------------------------------------------------------------

function addPoints(kmlDoc, lat, lon, options)
% Add points to the KML document.

% Verify the coordinates.
[lat, lon] = verifyCoordinates(lat, lon);

% Append the points to the KML document.
appendPoint(kmlDoc, lat, lon, options);

%--------------------------------------------------------------------------

function addStruct(kmlDoc, S, options, userSupplied)
% Add a geostruct to the KML document.

% Verify the geostruct.
position = 2;
S = checkgeostruct(S, position, mfilename);

% Currently, only Point or MultiPoint geometry is supported.
assert(~isempty(strmatch(lower(S(1).Geometry), {'point','multipoint'})), ...
    'map:kmlwrite:invalidGeometry', ...
    ['Function KMLWRITE expected a geostruct with ''Point'' or ''MultiPoint'' geometry\n', ...
    'but a geometry of ''%s'' was found instead.'], S(1).Geometry);
       
% The coordinates of the geostruct have not been verified. 
% The coordinates may contain single or multi-point data. 
% Validate the coordinates outside of the loop for speed.
lat = extractfield(S,'Lat');
lon = extractfield(S,'Lon');
verifyCoordinates(lat, lon);

% Make the default HTML table if required.
options = makeDefaultTable(S, options, userSupplied);

% Via addMultiPointStruct, add the point data to the document, using
% arrayfun to process each element in the geostruct. 
fcn = @(x, y)appendMultiPointStruct(kmlDoc, x, y);
arrayfun(fcn, S, reshape(options,size(S)));

%--------------------------------------------------------------------------

function options = makeDefaultTable(S, options, userSupplied)
% Create a default HTML table, if needed. The options Description field
% will contain the table.

% Obtain the attribSpec either from the options structure if supplied by
% the user or create a default one.
if ~isempty(options(1).Spec)
    attribSpec = options(1).Spec;   
else
    attribSpec = makeattribspec(S);
end

% A table needs to be created if the user did not supply a description or
% an attribute spec is supplied.  If an attribute spec is supplied, then
% userSupplied.Description is previously set to false.
needTable = ~userSupplied.Description;

if needTable
    % Convert the geostruct to a html cell array. The table is the
    % description field of a KML Placemark element which is located at the
    % specified coordinates.
    html = geostructToHTML(S, attribSpec);
    for k=1:numel(S)
        % Copy html to the description.
        options(k).Description = html{k};
    end
end
