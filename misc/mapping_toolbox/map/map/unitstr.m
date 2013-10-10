function str = unitstr(str,typeOfUnits)
%UNITSTR  Check unit strings or abbreviations
%
%   UNITSTR is obsolete and will be removed in a future release.  The
%   syntax STR = UNITSTR(STR,'times') has already been removed.
%
%   UNITSTR, with no arguments, displays a list of strings and
%   abbreviations, recognized by certain Mapping Toolbox functions,
%   for units of angle and length/distance.
%   
%   STR = UNITSTR(STR,'angles') checks for valid angle unit strings or
%   abbreviations.  If a valid string or abbreviation is found, it
%   is converted to a standardized, preset string.
%
%   STR = UNITSTR(STR,'distances') checks for valid length unit strings
%   or abbreviations.  If a valid string or abbreviation is found, it is
%   converted to a standardized, preset string.  Note that 'miles' and
%   'mi' are converted to 'statutemiles';  there is no way to specify
%   international miles in the UNITSTR function.

% Copyright 1996-2009 The MathWorks, Inc.
% $Revision: 1.12.4.4 $  $Date: 2009/05/14 17:05:35 $

warning('map:unitstr:obsolete', ...
    'Function %s is obsolete and will be removed in a future release of the toolbox.', ...
    upper(mfilename))

if nargin == 0
    unitstra
    unitstrd
else
    error(nargchk(2, 2, nargin, 'struct'))

    validTypesOfUnits = {'angles','distances','times'};
    k = find(strncmpi(typeOfUnits, validTypesOfUnits, numel(typeOfUnits)));
    if numel(k) == 1
        typeOfUnits = validTypesOfUnits{k};
    else
        error([getcomp ':' mfilename ':unsupportTypeOfUnits'], ...
            'Unrecognized type of units: %s', typeOfUnits);
    end

    switch typeOfUnits
        case 'angles'
            str = checkangleunits(str);

        case 'distances'
            str = unitstrd(str);

        case 'times'
            error([getcomp ':' mfilename ':obsoleteTimeUnits'], ...
                'UNITSTR no longer supports time units.')
    end
end

%-----------------------------------------------------------------------

function unitstra
% Display a list of supported angle unit strings

disp(' ');    disp('Recognized Angle Unit Strings')
disp(' ');    disp(['degrees'; 'radians'])
disp(' ');    disp('Recognized Angle Unit Abbreviations')
disp(' ');    disp(['deg   for degrees'; 'rad   for radians'])
