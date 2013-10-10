function mstruct = setOrigin(mstruct, origin)
% Set/reset the Origin property of a map projection structure.

% Copyright 2008 The MathWorks, Inc.
% $Revision: 1.1.6.1 $  $Date: 2008/05/14 22:01:34 $

% Validate origin
if ischar(origin) || numel(origin) > 3
    error(['map:' mfilename ':mapdispError'], ...
        'Incorrect Origin property')
end

if numel(origin) == 3 && ~isempty(mstruct.fixedorient)
    warning(['map:', mfilename, ':FixedOrientationWarn'], ...
        'Fixed orientation supplied for this projection.')
    mstruct.origin = origin;
elseif strcmp(mstruct.mapprojection, 'utm')
    warning(['map:', mfilename, ':UTMZoneWarn'], ...
        ['Ignoring origin value. ' ...
        'The origin of a UTM system is determined by its Zone property.'])
elseif strcmp(mstruct.mapprojection, 'ups')
    warning(['map:', mfilename, ':UPSZoneWarn'], ...
        ['Ignoring origin value. ' ...
        'The origin of a UPS system is determined by its Zone property.'])
elseif isscalar(origin)
    % Interpret scalar value as origin longitude.
    mstruct.origin = [NaN origin];
else
    mstruct.origin = origin;
end
