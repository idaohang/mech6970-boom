function tf = okToSetMapLimits(mstruct,propname)
% True unless the projection type and/or orientation angle precludes
% controlling frame limits of a map projection structure by specifying
% MapLatLimit and MapLonLimit vectors.

% Copyright 2008-2009 The MathWorks, Inc.
% $Revision: 1.1.6.2 $  $Date: 2009/05/14 17:06:39 $

tf = true;
if mprojIsExceptional(mstruct.mapprojection)
    warning(['map:', mfilename, ':ignoring' propname '1'],...
        ['Ignoring value of %s. %s cannot be set ', ...
        'for the %s projection.'], propname, propname, mstruct.mapprojection)
    tf = false;
elseif ~isempty(mstruct.origin)
    origin = mstruct.origin;
    nonzeroOrientationAngle = (numel(origin) >= 3) && (origin(3) ~= 0);    
    if nonzeroOrientationAngle
        warning(['map:', mfilename, ':ignoring' propname '2'],...
            'Ignoring value of %s due to nonzero orientation angle.', ...
            propname)
        tf = false;
    elseif ~mprojIsAzimuthal(mstruct.mapprojection)
        originlat = mstruct.origin(1);
        nonzeroOriginLat = ~isnan(originlat) && (originlat ~= 0);
        if nonzeroOriginLat && mprojIsRotational(mstruct.mapprojection)
            warning(['map:', mfilename, ':ignoring' propname '3'],...
                ['Ignoring value of %s due to use of nonzero origin ', ...
                'latitude with the %s projection.'], ...
                propname, mstruct.mapprojection)
            tf = false;
        end
    end
end

%-----------------------------------------------------------------------

function tf = mprojIsExceptional(mapprojection)

exceptions = {...
    'globe', ...    % No need for map limits -- always covers entire planet
    'cassini', ...  % Always in a transverse aspect
    'wetch', ...    % Always in a transverse aspect
    'bries'};       % Always in an oblique aspect
tf = any(strcmpi(mapprojection,exceptions));

%-----------------------------------------------------------------------

function tf = mprojIsRotational(mapprojection)

tf = ~any(strcmp(mapprojection, {...
       'utm', 'tranmerc', 'lambertstd', 'cassinistd', ...
       'eqaconicstd', 'eqdconicstd', 'polyconstd'}));
