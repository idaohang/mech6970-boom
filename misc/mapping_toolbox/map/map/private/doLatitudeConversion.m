function latout = doLatitudeConversion(function_name,from,to,varargin)
% DOLATITUDECONVERSION  Compute engine for old latitude conversion functions.
%
%   See also:  AUT2GEOD, CEN2GEOD, CNF2GEOD, ISO2GEOD, GEOD2AUT,
%              GEOD2CEN, GEOD2CNF, GEOD2ISO, GEOD2PAR, GEOD2REC,
%              PAR2GEOD, REC2GEOD.

% Copyright 1996-2010 The MathWorks, Inc.
% $Revision: 1.1.6.4 $  $Date: 2010/09/24 14:33:50 $

try
    error(nargchk(1,3,numel(varargin),'struct'))
catch exception
    exception.throwAsCaller
end

[latin,ellipsoid,units] = parseLatConvInputs(function_name,varargin{:});
latin = toRadians(units, latin);
latout = convertlat(ellipsoid,latin,from,to,'nocheck');
latout = fromRadians(units, latout);

%--------------------------------------------------------------------------

function [latin,ellipsoid,units] = parseLatConvInputs(function_name,varargin)

latin = varargin{1};
switch(numel(varargin))
    case 1
        ellipsoid = almanac('earth','geoid');
	    units = 'degrees';
    case 2
        if ischar(varargin{2});
            ellipsoid = almanac('earth','geoid');
            units = varargin{2};
        else
            ellipsoid = varargin{2};
            units = 'degrees';
        end
    case 3
        ellipsoid = varargin{2};
        units = varargin{3};
end

ellipsoid = checkellipsoid(ellipsoid,function_name,'ELLIPSOID',2);
