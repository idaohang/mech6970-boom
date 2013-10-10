function html = geostructToHTML(S, attribSpec)
%GEOSTRUCTTOHTML Convert geostruct to HTML
%
%   html = geostructToHTML(S, ATTRIBSPEC) outputs an HTML table to the
%   string cell array HTML for each entry in the geostruct S by applying
%   the attribute specification, attribSpec.

% Copyright 2007 The MathWorks, Inc.
% $Revision: 1.1.6.3 $  $Date: 2007/07/31 20:03:15 $

% Validate the attribute spec.
attribSpec = validateAttributeSpec(attribSpec, S);

% Obtain the field names of the attribute structure to be used in the
% table.
rowHeaderNames = getAttributeLabelsFromSpec(attribSpec);

% Convert the structure to a string cell array.
c = convertStructToStringCell(S, attribSpec);

% Convert each element of the geostruct to a HTML table. 
html = cell(1,numel(S));
for k=1:numel(S)
    % Convert the cell array to an HTML table.
    html{k} = makeTable(rowHeaderNames, c(:,k));
end

%--------------------------------------------------------------------------

function attributeLabels = getAttributeLabelsFromSpec(attribSpec)
% Obtain the table names from the attribute spec.  attributeLabels is a
% cell array containing the string names for each attribute.

% Assign the field names of the desired attributes.
specFieldNames = fieldnames(attribSpec);

% Assign the number of desired attributes.
numSpecFields = numel(specFieldNames);

% Allocate space for the attribute labels.
attributeLabels = cell(1,numSpecFields);

% Obtain the names from the spec.
for m=1:numSpecFields
    specFieldName = specFieldNames{m};
    attributeLabels{m} = attribSpec.(specFieldName).AttributeLabel;   
end

%--------------------------------------------------------------------------

function c = convertStructToStringCell(S, attribSpec)
% Convert the structure array S to a string cell array by applying the 
% format in the structure attribSpec.

% Apply the spec to the structure S.
A = getAttributeStruct(S, attribSpec);
A = applyAttributeSpec(A, attribSpec);

% Convert the structure array to a cell array.
c = struct2cell(A);

%--------------------------------------------------------------------------

function A = getAttributeStruct(S, spec)
% Get the attributes from the geostruct S. A is a struct the same size as
% S containing a subset of the fields of S -- all valid attribute fields.
% The order of the fields of A matches the order in spec.

% Obtain the fieldnames of S.
sFieldNames = fieldnames(S);

% Obtain an index of field names that are not in the attribute spec, but
% are contained in S. These fields will be removed to form an attribute
% structure.
nonSpecIndex = ~ismember(sFieldNames, fieldnames(spec));

% Remove all non-attribute fields in S and all fields that are not in the
% attribute spec.
A = rmfield(S, sFieldNames(nonSpecIndex));

% A contains field names in the same order as spec.
A = orderfields(A, spec);

%--------------------------------------------------------------------------

function A = applyAttributeSpec(A, attribSpec)
% Apply the attribute spec, attribSpec, to the attribute structure A.

% Assign the field names of the desired attributes.
fieldNames = fieldnames(attribSpec);

% Loop through each element of the attribute structure and convert each
% field value to a string.
A = arrayfun(@(a) attributeToString(a, fieldNames, attribSpec), A);

%--------------------------------------------------------------------------

function A = attributeToString(A,  fieldNames, attribSpec)
% Convert field values in attribute structure to string values.

% For each field of A, convert the value of the field to a string, using
% the format from the attribute spec.
for k=1:numel(fieldNames)
    fieldName = fieldNames{k};
    A.(fieldName) = num2str(A.(fieldName), attribSpec.(fieldName).Format);
end

%--------------------------------------------------------------------------

function attribSpec = validateAttributeSpec(attribSpec, S)

% If attribSpec is empty, make sure it's an empty struct.
if isempty(attribSpec)
    attribSpec = struct([]);
    % No need to check anything else, including the attribute fields of S.
    return  
end

% Make sure that attribSpec is a structure.
assert(isstruct(attribSpec), ...
    'map:geostructToHTML:attribSpecNotAStructure', ...
    'Function KMLWRITE expected the attribute spec to be a structure.');

% Make sure that attribSpec is a scalar (or empty).
assert(numel(attribSpec) == 1, ...
    'map:geostructToHTML:nonScalarAttribSpec', ...
    'Function KMLWRITE expected the attribute spec to be a scalar (or empty) structure.');

% Validate S, then make sure that attribSpec and S are mutually consistent.
defaultspec = makeattribspec(S);  % Validates attribute values in S
attributeNamesInS = fields(defaultspec);
attributeNamesInattribSpec = fields(attribSpec);

% Check 1:  Every attribute in attribSpec must be present in S.
missingAttributes = setdiff(attributeNamesInattribSpec,attributeNamesInS);
assert(isempty(missingAttributes), ...
    'map:geostructToHTML:missingAttributes', ...
    'The attribute spec specifies attributes that are missing from the geostruct S.');
 
% Check 2:  Field types in attribSpec must be consistent with S.
% While in this loop, use the default to fill in any fields missing from 
% attribSpec.
for k = 1:numel(attributeNamesInattribSpec)
    attributeName = attributeNamesInattribSpec{k};
    fSpec = attribSpec.(attributeName);
    fDefault = defaultspec.(attributeName);
    
    if ~isfield(fSpec,'AttributeLabel')
        attribSpec.(attributeName).AttributeLabel = fDefault.AttributeLabel;
    end

    if ~isfield(fSpec,'Format')
        attribSpec.(attributeName).Format = fDefault.Format;
    end 
end

%--------------------------------------------------------------------------

function  html = makeTable(names, values)
% Create a cell array containing HTML embedded tags representing a table.
% The table is two-column with name, value pairs. 

% NAMES is a string cell array containing the names for the first column in
% the table. VALUES is a string cell array containing the values for the
% second column in the table.  HTML is a char array containing embedded
% HTML tags defining a table. 

if ~isempty(names)
    numRows = numel(names);
    html = cell(numRows+2,1);

    html{1} = sprintf('<html><table border="1">\n');
    rowFmt = '<tr><td>%s</td><td>%s</td></tr>\n';
    for k = 1:numRows
        html{k+1} = sprintf(rowFmt, names{k}, values{k});
    end
    html{numRows+2} = sprintf('</table><br></br></html>\n');
    html = char([html{:}]);
else
    html = ' ';
end
