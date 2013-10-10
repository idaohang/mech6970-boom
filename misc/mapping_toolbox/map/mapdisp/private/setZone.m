function mstruct = setZone(mstruct, zone)
% Update a map projection structure given a new value for its Zone property.

% Copyright 2008 The MathWorks, Inc.
% $Revision: 1.1.6.1 $  $Date: 2008/05/14 22:01:35 $

switch mstruct.mapprojection
    case 'ups'
        if ~any(strcmp(zone,{'north','south'}))
            error(['map:' mfilename ':mapdispError'], ...
                'Incorrect %s property. Recognized zones are %s and %s.', ...
                '''zone''','''north''','''south''')
        end
        mstruct.zone = zone;
        mstruct = feval(mstruct.mapprojection,mstruct);
        mstruct.maplatlimit = [];
        mstruct.flatlimit = [];
        mstruct.origin = [];
        mstruct.mlabelparallel = [];
        mstruct.mlinelimit = [];
    case 'utm'
        [latlim,lonlim,txtmsg] = utmzone(zone);
        if ~isempty(txtmsg)
            error(['map:' mfilename ':mapdispError'], ...
                ['Incorrect %s property. Recognized zones are integers\n' ...
                'from 1 to 60 or numbers followed by letters from C to X.'], ...
                '''zone''')
        else
            mstruct.zone = upper(zone);
            mstruct.flatlimit = fromDegrees(mstruct.angleunits,latlim);
            mstruct.flonlimit = fromDegrees(mstruct.angleunits,[-3 3]);
            mstruct.origin = fromDegrees(mstruct.angleunits,[0 min(lonlim)+3 0]);
            mstruct.maplatlimit = [];
            mstruct.maplonlimit = [];

            mstruct.mlinelocation = [];
            mstruct.plinelocation = [];
            mstruct.mlabellocation = [];
            mstruct.plabellocation = [];
            mstruct.mlabelparallel = [];
            mstruct.plabelmeridian = [];
            mstruct.falsenorthing = [];
        end
    otherwise
        error(['map:' mfilename ':mapdispError'], ...
            'ZONE cannot be specified for this projection')
end
