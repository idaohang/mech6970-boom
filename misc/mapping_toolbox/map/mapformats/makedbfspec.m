function dbfspec = makedbfspec(S)
%MAKEDBFSPEC DBF specification from geographic data structure
%  
%   DBFSPEC = MAKEDBFSPEC(S) analyzes a geographic data structure array,
%   S (a geostruct or mapstruct), and constructs a DBF specification
%   suitable for use with SHAPEWRITE.  You can modify DBFSPEC, then pass
%   it to SHAPEWRITE to exert control over which attribute fields are
%   written to the DBF component of the shapefile, the field-widths, and
%   the precision used for numerical values.
%
%   DBFSPEC is a scalar MATLAB structure with two levels.  The top level
%   consists of a field for each attribute in S.  Each of these fields,
%   in turn, contains a scalar structure with a fixed set of four fields:
%
%   FieldName          The field name to be used within the DBF file.  This
%                      will be identical to the name of the corresponding
%                      attribute, but may modified prior to calling
%                      SHAPEWRITE.  This might be necessary, for example,
%                      because you want to use spaces your DBF field names,
%                      but the attribute fieldnames in S must be valid
%                      MATLAB variable names and cannot have spaces
%                      themselves.
%
%   FieldType          The field type to be used in the file, either 'N'
%                      (numeric) or 'C' (character).
%
%   FieldLength        The number of bytes that each instance of the field
%                      will occupy in the file.
%
%   FieldDecimalCount  The number of digits to the right of the decimal
%                      place that are kept in a numeric field. Zero for
%                      integer-valued fields and character fields. The
%                      default value for non-integer numeric fields is 6.
%
%   Example
%   -------
%   % Import a shapefile representing a small network of road segments,
%   % and construct a DBF specification.
%   s = shaperead('concord_roads')
%   dbfspec = makedbfspec(s)
%
%   % Modify the DBF spec to (a) eliminate the 'ADMIN_TYPE' attribute, (b)
%   % rename the 'STREETNAME' field to 'Street Name', and (c) reduce the
%   % number of decimal places used to store road lengths.
%   dbfspec = rmfield(dbfspec,'ADMIN_TYPE')
%   dbfspec.STREETNAME.FieldName = 'Street Name';
%   dbfspec.LENGTH.FieldDecimalCount = 1;
%
%   % Export the road network back to a modified shapefile.  (Actually,
%   % only the DBF component will be different.)
%   shapewrite(s, 'concord_roads_modified', 'DbfSpec', dbfspec)
%
%   % Verify the changes you made.  Notice the appearance of
%   % 'Street Name' in the field names reported by SHAPEINFO, the absence
%   %  of the 'ADMIN_TYPE' field, and the reduction in the precision of the
%   %  road lengths.
%   info = shapeinfo('concord_roads_modified')
%   {info.Attributes.Name}
%   r = shaperead('concord_roads_modified')
%   s(33).LENGTH
%   r(33).LENGTH
%
%   See also: SHAPEINFO, SHAPEWRITE.

% Copyright 2003-2009 The MathWorks, Inc.  
% $Revision: 1.1.6.11 $  $Date: 2009/09/03 05:17:51 $

% Make sure S is non-empty.
if isempty(S)
    eid = sprintf('%s:%s:emptyGeostruct', getcomp, mfilename);
    error(eid, 'S must be a non-empty geographic data structure array.')
end

% Make sure that S is a struct array.
if ~isstruct(S)
    eid = sprintf('%s:%s:notAStruct', getcomp, mfilename);
    error(eid, 'S must be a structure array.')
end

% Determine what types of fields to write for each attribute.
attributeNames = fields(S);
[~,fIndex] = ...
    setxor(attributeNames, {'Geometry','BoundingBox','X','Y','Lat','Lon'});
attributeNames = attributeNames(sort(fIndex));

% Default to six digits to the right of the decimal point.
defaultDecimalPrecision = 6;

% In this version we support only types 'N' and 'C'.
for k = 1:numel(attributeNames)
    dataClass = class(S(1).(attributeNames{k}));
    switch(dataClass)
        
        case 'double'
            fieldType = 'N';
            try
                v = [S.(attributeNames{k})];
            catch exception
                % Issue our own error, unless MATLAB has run out of memory.
                if strcmp(exception.identifier,'MATLAB:nomem')
                    rethrow(exception)
                else
                    error([getcomp ':' mfilename ':inconsistentAttributeData1'], ...
                        'Inconsistent data in attribute field: %s.', ...
                        attributeNames{k})
                end
            end
            
            if ~isa(v,'double')
                error([getcomp ':' mfilename ':nondoubleNumericAttribute'], ...
                    ['Attribute field %s of geographic data structure S', ...
                    ' contains at least one non-double value.'], ...
                    attributeNames{k})
            end
            
            if numel(v) ~= numel(S)
                error([getcomp ':' mfilename ':nonscalarAttributeValue'],...
                    ['Attribute field %s of geographic data structure S', ...
                    ' contains at least one value that is not a scalar double.'], ...
                    attributeNames{k})
            end
            
            if any(isinf(v)) || any(~isreal(v))
                error([getcomp ':' mfilename ':attributeNotFiniteReal'],...
                    ['Numerical attributes of geographic data structure S', ...
                    ' must be finite and real.'])
            end
                
            if all(v == 0)
                numRightOfDecimal = 0;
                fieldLength = 2;
            else
                numLeftOfDecimal = max(1, 1 + floor(log10(max(abs(v)))));
                if all(v == floor(v))
                    numRightOfDecimal = 0;
                    fieldLength = 1 + numLeftOfDecimal;
                else
                    numRightOfDecimal = defaultDecimalPrecision;
                    fieldLength = 1 + numLeftOfDecimal + 1 + numRightOfDecimal;
                end
            end
            minNumericalFieldLength = 3;  % Large enough to hold 'NaN'
            fieldLength = max(fieldLength, minNumericalFieldLength);

        case 'char'
            fieldType = 'C';
            numRightOfDecimal = 0;
            
            % Obtain the attribute values from the structure array.
            values = {S.(attributeNames{k})};
              
            % Calculate the required field length for the attribute.
            maxFieldLength = ...
                calculateMaxFieldLength(attributeNames{k}, values);
            minCharFieldLength = 2;
            fieldLength = max(maxFieldLength, minCharFieldLength);
                       
        otherwise
            wid = sprintf('%s:%s:unsupportedDataClass', getcomp, mfilename);
            warning(wid, 'Omitting unsupported data class: %s', dataClass);
            continue;
    end    
    
    dbfspec.(attributeNames{k}) = struct(...
        'FieldName', attributeNames{k},...
        'FieldType', fieldType,...
        'FieldLength', fieldLength,...
        'FieldDecimalCount', numRightOfDecimal);
end

% Handle both:
%    (1) a geostruct or mapstruct without an attributes and
%    (2) a geostruct or mapstruct in which every attribute contains
%        an unsupported data class.
if ~exist('dbfspec','var')
    dbfspec = struct([]);
end

%--------------------------------------------------------------------------

function maxFieldLength = calculateMaxFieldLength(attributeName, values)
% Calculate a safe upper bound on the number of bytes per field required to
% contain the characters in the cell array, VALUES.  attributeName is the
% name of the attribute.

% Validate that the cell array, VALUES, contains only char characters.
try
    charValues = char(values);
catch exception
    % Issue our own error.
    if strcmp(exception.identifier,'MATLAB:char:CellsMustContainChars')
        error([getcomp ':' mfilename ':inconsistentAttributeData2'], ...
            ['Attribute field %s of geographic data structure %s', ...
            ' contains at least one value that is not a character string.'],...
            attributeName, 'S')      
    else
        rethrow(exception)
    end
end
            
% Determine a safe upper bound for the number of bytes per field required
% to contain the number of characters. For efficiency, use charValues to
% determine if the input, VALUES, is a cell array containing all ASCII
% characters. 
if isASCII(charValues)
    % charValues contains all ASCII characters. A safe upper bound is the
    % size of the rows.
    maxFieldLength = size(charValues, 2);
else
    % The array contains non-ASCII unicode characters. Calculate
    % maxFieldLength as the maximum number of bytes required to hold any
    % element in VALUES. For each element of the cell array, convert the
    % element to native representation and count the resultant number of
    % bytes. Set maxFieldLength as the maximum number of bytes of any
    % element in the cell array.
    numNativeBytesFcn = @(x)numel(unicode2native(x));
    numelCell = cellfun(numNativeBytesFcn, values, 'UniformOutput', false);
    maxFieldLength = max(cell2mat(numelCell));
end

%--------------------------------------------------------------------------

function tf = isASCII(c)
% Return true if the character array, C, contains all ASCII characters.

tf = isequal(char(uint8(c)), c);
