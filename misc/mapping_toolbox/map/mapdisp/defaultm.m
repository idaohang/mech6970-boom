function [mstruct,msg] = defaultm(varargin)
%DEFAULTM Initialize or reset map projection structure
%
%   MSTRUCT = DEFAULTM(PROJECTION) initializes a map projection
%   structure.  PROJECTION is a string containing the name of a
%   Mapping Toolbox projection function.
%
%   MSTRUCT = DEFAULTM(MSTRUCT) checks an existing map projection
%   structure, sets empty properties, and adjusts dependent properties.
%   The Origin, FLatLimit, FLonLimit, MapLatLimit, and MapLonLimit
%   properties may be adjusted for compatibility with each other and
%   with the MapProjection and (in the case of UTM or UPS) Zone
%   properties.
%
%   See also AXESM, MAPLIST, MAPS, MFWDTRAN, MINVTRAN, SETM.

% Copyright 1996-2010 The MathWorks, Inc.
% $Revision: 1.9.4.4 $  $Date: 2010/03/04 16:22:37 $

% Obsolete Syntax
% ---------------
%  [MSTRUCT, MSG] = DEFAULTM(...) returns a string indicating any error
%  encountered.
if nargout > 1
    warnObsoleteMSGSyntax(mfilename)
    msg = '';
end

% Parse inputs, then either initialize or reset.
if nargin == 0 || ischar(varargin{1})
    if nargin == 0
        % MSTRUCT = DEFAULTM
        mapprojection = [];
    else
        % MSTRUCT = DEFAULTM(PROJECTION)
        mapprojection = varargin{1};
    end
    mstruct = initmstruct;
    mstruct.mapprojection = mapprojection;
    if ~isempty(mstruct.mapprojection)
        mstruct = feval(mstruct.mapprojection,mstruct);
    end
else
    % MSTRUCT = DEFAULTM(MSTRUCT)
    mstruct = varargin{1};
        
    % Pre-process key, interrelated properties in the same order as in
    % AXESM.  Start by saving their input values in local variables.
    zone          = mstruct.zone;
    origin        = mstruct.origin;
    maplatlimit   = mstruct.maplatlimit;
    maplonlimit   = mstruct.maplonlimit;
        
    if ~isempty(zone)
        mstruct = setZone(mstruct, zone);
    end
    
    if ~isempty(origin)
        mstruct = setOrigin(mstruct, origin);
    end

    if ~isempty(maplatlimit)
        mstruct = setMapLatLimit(mstruct, maplatlimit);
    end

    if ~isempty(maplonlimit)
        mstruct = setMapLonLimit(mstruct, maplonlimit);
    end

    %  Remove possible NaN left by setOrigin or setMapLatLimit.
    mstruct.origin(isnan(mstruct.origin)) = 0;

    % Ensure that all properties are assigned reasonable, consistent
    % values, following the same procedure as in AXESM.
    mstruct = resetmstruct(mstruct);
end
