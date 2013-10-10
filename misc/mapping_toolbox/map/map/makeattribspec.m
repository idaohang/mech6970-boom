function attribspec = makeattribspec(S)
%MAKEATTRIBSPEC Attribute specification from geographic data structure
%
%   ATTRIBSPEC = MAKEATTRIBSPEC(S) analyzes a geographic data structure S
%   and constructs an attribute specification suitable for use with
%   KMLWRITE.  KMLWRITE, given geostruct input, constructs a HTML table
%   that consists of a label for the attribute in the first column and the
%   string value of the attribute in the second column.  You can modify
%   ATTRIBSPEC, then pass it to KMLWRITE to exert control over which
%   geostruct attribute fields are written to the HTML table and the format
%   of the string conversion.
%
%   ATTRIBSPEC is a scalar MATLAB structure with two levels.  The top level
%   consists of a field for each attribute in S.  Each of these fields, in
%   turn, contains a scalar structure with a fixed pair of fields:
%
%   AttributeLabel     A string that corresponds to the name of the
%                      attribute field in the geographic data structure.
%                      With KMLWRITE, the string is used to label the
%                      attribute in the first column of the HTML table.
%                      The string may be modified prior to calling
%                      KMLWRITE.  You might modify an attribute label,
%                      for example, because you want to use spaces in
%                      your HTML table, but the attribute fieldnames in
%                      S must be valid MATLAB variable names and cannot
%                      have spaces themselves.
%
%   Format             The sprintf format character string that converts 
%                      the attribute value to a string.
%
%   Example
%   -------
%   % Import a shapefile representing tsunami (tidal wave) events reported 
%   % over several decades, tagged geographically by source location.
%   s = shaperead('tsunamis', 'UseGeoCoords', true);
%
%   % Construct an attribute specification.
%   attribspec = makeattribspec(s);
%
%   % Modify the attribute spec to:
%   % (a) Display Max_Height, Cause, Year, Location, and Contry attributes 
%   % (b) Rename the 'Max_Height' field to 'Maximum Height' 
%   % (c) Highlight each attribute label with a bold font 
%   % (d) Set to zero the number of decimal places used to display Year
%   % (e) We have independent knowledge that the height units are meters, 
%   %     so we will add that to the Height format specifier
%
%   desiredAttributes = ...
%      {'Max_Height', 'Cause', 'Year', 'Location', 'Country'};
%   allAttributes = fieldnames(attribspec);
%   attributes = setdiff(allAttributes, desiredAttributes);
%   attribspec = rmfield(attribspec, attributes);
%   attribspec.Max_Height.AttributeLabel = '<b>Maximum Height</b>';
%   attribspec.Max_Height.Format = '%.1f Meters';
%   attribspec.Cause.AttributeLabel = '<b>Cause</b>';
%   attribspec.Year.AttributeLabel = '<b>Year</b>';
%   attribspec.Year.Format = '%.0f';
%   attribspec.Location.AttributeLabel = '<b>Location</b>';
%   attribspec.Country.AttributeLabel = '<b>Country</b>';
% 
%   % Export the selected attributes and source locations to a KML file. 
%   filename = 'tsunami.kml';
%   kmlwrite(filename, s, 'Description', attribspec, 'Name', {s.Location})
%
%   See also KMLWRITE, MAKEDBFSPEC, SHAPEWRITE.

% Copyright 2007 The MathWorks, Inc.
% $Revision: 1.1.6.5 $  $Date: 2008/01/15 18:53:54 $

% Make sure S is non-empty.
assert(~isempty(S), ...
    'map:makeattribspec:emptyGeostruct', ...
    'S must be a non-empty geostruct array.');

% Make sure that S is a struct array.
assert(isstruct(S), ...
    'map:makeattribspec:notAStruct', ...
    'Geostruct S must be a structure array.');

% Determine attribute fields.
attributeNames = fields(S);
[t1,fIndex] = setxor(attributeNames, ...
    {'Geometry','BoundingBox','X','Y','Lat','Lon'});
attributeNames = attributeNames(sort(fIndex));

for k = 1:numel(attributeNames)
    v = extractfield(S, attributeNames{k});

    if  ~iscell(v) 
        % Validate numeric attributes.
        assert(numel(v) == numel(S), ...
            'map:makeattribspec:nonscalarAttributeValue', ...
            ['Attribute field %s of geostruct S\n', ...
            'contains at least one value that is not a scalar.'],...
            attributeNames{k});

        assert(all(~isinf(v)) && all(isreal(v)), ...
            'map:makeattribspec:attributeNotFiniteReal', ...
            'Numerical attributes of geostruct S must be finite and real.');

        format = '%.15g';
        attribspec.(attributeNames{k}) = struct(...
            'AttributeLabel', attributeNames{k},...
            'Format', format);
    else
        dataClass = unique(cellfun(@class, v, 'UniformOutput', false));
        if numel(dataClass) == 1 && ~isequal(dataClass{1}, 'char')
            warning('map:makeattribspec:unsupportedDataClass', ...
                'Omitting unsupported data class: %s', dataClass{1});
        else
            mustBeChar = (numel(dataClass) == 1) && ...
                isequal(dataClass{1}, 'char');
            assert(mustBeChar, ...
                'map:makeattribspec:inconsistentDataClass', ...
                'Inconsistent data classes in attribute field: %s', ...
                attributeNames{k});

            format = '%s';
            attribspec.(attributeNames{k}) = struct(...
                'AttributeLabel', attributeNames{k},...
                'Format', format);
        end
    end
end

% Return empty if there are no attributes or attribspec is unassigned.
if isempty(attributeNames) || ~exist('attribspec','var')
    attribspec = struct([]);
end
