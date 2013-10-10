function code = projcode(id)
%PROJCODE Obtain a projection code structure
%
%   CODE = PROJCODE(ID) returns a projection code given an ID. ID may be a
%   string or scalar.  ID may be either an mstruct mapprojection name or a
%   GeoTIFF projection name if ID begins with 'CT_'.  The output CODE
%   structure contains the following fields:
%      index:         sequential index number
%      mapprojection: mstruct mapprojection name 
%      projname:      PROJ4 projection name 
%      CTProjection:  GeoTIFF projection name  
%      CTcode:        GeoTIFF projection number        
%
%   If the ID cannot be determined, 'Unknown' will be set for all names and
%   the maximum index number will be returned.
%
%   Example
%   -------
%      % Obtain the code structure for the tranmerc projection
%      %  using all parameter types.
%      code = projcode('CT_TransverseMercator');
%      code = projcode('tranmerc');
%      code = projcode(1);
%
%   See also GEOTIFF2MSTRUCT, GEOTIFFINFO.

% Copyright 1996-2006 The MathWorks, Inc.
% $Revision: 1.1.6.4 $  $Date: 2006/06/15 20:13:45 $

% Obtain the conversion data
geopairs = GeoValuePairs;
prjpairs = ProjValuePairs;
mappairs = MapValuePairs;
prjnames = NameValuePairs;
projMax  = numel(MapValuePairs);

if (ischar(id)) 
   % Process the Projection name
   if (isempty(strmatch('CT_',id)))
      % mstruct name
      pairs = mappairs;
   else
      % GeoTIFF name
      pairs = geopairs;
   end

   % Find the name and number
   projnum = projMax;
   for k=1:projMax
      if (strcmp(pairs{k}.name, id))
         projnum = pairs{k}.number;
         break;
      end
   end
else
   % Input is a number
   projnum = id;
   try
      checkinput(id, {'numeric'}, ...
                {'real' 'scalar' 'positive', 'integer'}, ...
                mfilename, 'ID', 1);
   catch
      projnum = projMax;
   end
   if projnum > projMax
     projnum = projMax;
   end
end

% Set the projection code given the projection number
code.index = projnum;
code.mapprojection = mappairs{projnum}.name;
code.projname      = prjpairs{projnum}.name;
code.CTProjection  = geopairs{projnum}.name;
code.CTcode        = geopairs{projnum}.number;
code.Name          = prjnames{projnum}.name;

%--------------------------------------------------------------------------
function pair = ValuePair(name, number)
pair.name   = name;
pair.number = number;

%--------------------------------------------------------------------------
function pairs = MapValuePairs
% mstruct projection name pairs
%  unmatched pairs return Unknown
pairs = {ValuePair('tranmerc',1), ...
         ValuePair('Unknown', 2), ...
         ValuePair('Unknown', 3), ...
         ValuePair('Unknown', 4), ...
         ValuePair('Unknown', 5), ...
         ValuePair('Unknown', 6), ...   
         ValuePair('mercator',7), ...
         ValuePair('lambertstd', 8), ...
         ValuePair('lambertstd', 9), ...
         ValuePair('eqaazim', 10), ...
         ValuePair('eqaconicstd',11), ...
         ValuePair('eqdazim', 12), ...
         ValuePair('eqdconicstd',13), ...
         ValuePair('Unknown', 14), ...
         ValuePair('ups',     15), ...
         ValuePair('stereo',  16), ...   
         ValuePair('eqdcylin',17), ...
         ValuePair('cassinistd', 18), ...
         ValuePair('gnomonic',19), ...
         ValuePair('miller',  20), ...
         ValuePair('ortho',   21), ...
         ValuePair('polyconstd', 22), ...
         ValuePair('robinson',23), ...
         ValuePair('sinusoid',24), ...
         ValuePair('vgrint1', 25), ...
         ValuePair('Unknown', 26), ...
         ValuePair('Unknown', 27), ...
         ValuePair('Unknown', 32767)};

%--------------------------------------------------------------------------
function pairs = ProjValuePairs
% PROJ4 projection name pairs
%  unmatched pairs return Unknown
pairs = {ValuePair('tmerc',   1), ...
         ValuePair('Unknown', 2), ...
         ValuePair('omerc',   3), ...
         ValuePair('Unknown', 4), ...
         ValuePair('Unknown', 5), ...
         ValuePair('Unknown', 6), ...   
         ValuePair('merc',    7), ...
         ValuePair('lcc',     8), ...
         ValuePair('lcc',     9), ...
         ValuePair('laea',    10), ...
         ValuePair('aea',     11), ...
         ValuePair('aeqd',    12), ...
         ValuePair('eqdc',    13), ...
         ValuePair('stere',   14), ...
         ValuePair('stere',   15), ...
         ValuePair('stere',   16), ...   
         ValuePair('eqc',     17), ...
         ValuePair('cass',    18), ...
         ValuePair('gnom',    19), ...
         ValuePair('mill',    20), ...
         ValuePair('ortho',   21), ...
         ValuePair('poly',    22), ...
         ValuePair('robin',   23), ...
         ValuePair('sinu',    24), ...
         ValuePair('vandg',   25), ...
         ValuePair('Unknown', 26), ...
         ValuePair('Unknown', 27), ...
         ValuePair('Unknown', 32767)};

%--------------------------------------------------------------------------
function pairs = GeoValuePairs
% GeoTIFF projection name pairs
%  Unknown is 32767
pairs = {ValuePair('CT_TransverseMercator',1), ...
         ValuePair('CT_TransvMercator_Modified_Alaska',2), ...
         ValuePair('CT_ObliqueMercator',3), ...
         ValuePair('CT_ObliqueMercator_Laborde',4), ...
         ValuePair('CT_ObliqueMercator_Rosenmund',5), ...
         ValuePair('CT_ObliqueMercator_Spherical',6), ...   
         ValuePair('CT_Mercator',7), ...
         ValuePair('CT_LambertConfConic_2SP',8), ...
         ValuePair('CT_LambertConfConic_1SP',9), ...
         ValuePair('CT_LambertAzimEqualArea',10), ...
         ValuePair('CT_AlbersEqualArea',11), ...
         ValuePair('CT_AzimuthalEquidistant',12), ...
         ValuePair('CT_EquidistantConic',13), ...
         ValuePair('CT_Stereographic',14), ...
         ValuePair('CT_PolarStereographic',15), ...
         ValuePair('CT_ObliqueStereographic',16), ...   
         ValuePair('CT_Equirectangular',17), ...
         ValuePair('CT_CassiniSoldner',18), ...
         ValuePair('CT_Gnomonic',19), ...
         ValuePair('CT_MillerCylindrical',20), ...
         ValuePair('CT_Orthographic',21), ...
         ValuePair('CT_Polyconic',22), ...
         ValuePair('CT_Robinson',23), ...
         ValuePair('CT_Sinusoidal',24), ...
         ValuePair('CT_VanDerGrinten',25), ...
         ValuePair('CT_NewZealandMapGrid',26), ...
         ValuePair('CT_TransvMercator_SouthOrientated',27), ...
         ValuePair('Unknown',32767)};

%--------------------------------------------------------------------------
function pairs = NameValuePairs
% GeoTIFF projection name pairs
%  Unknown is 32767
pairs = {ValuePair('Transverse Mercator',1), ...
         ValuePair('Transverse Mercator (Modified Alaska)',2), ...
         ValuePair('Oblique Mercator',3), ...
         ValuePair('Laborde Oblique Mercator',4), ...
         ValuePair('Rosenmund Oblique Mercator',5), ...
         ValuePair('Oblique Mercator Spherical',6), ...   
         ValuePair('Mercator',7), ...
         ValuePair('Lambert Conformal Conic (2SP)',8), ...
         ValuePair('Lambert Conformal Conic (1SP)',9), ...
         ValuePair('Lambert Azimuthal Equal Area',10), ...
         ValuePair('Albers Equal-Area Conic',11), ...
         ValuePair('Azimuthal Equidistant',12), ...
         ValuePair('Equidistant Conic',13), ...
         ValuePair('Stereographic',14), ...
         ValuePair('Polar Stereographic',15), ...
         ValuePair('Oblique Stereographic',16), ...   
         ValuePair('Equirectangular',17), ...
         ValuePair('Cassini-Soldner',18), ...
         ValuePair('Gnomonic',19), ...
         ValuePair('Miller Cylindrical',20), ...
         ValuePair('Orthographic',21), ...
         ValuePair('Polyconic',22), ...
         ValuePair('Robinson',23), ...
         ValuePair('Sinusoidal',24), ...
         ValuePair('Van der Grinten',25), ...
         ValuePair('New Zealand Map Grid',26), ...
         ValuePair('Transverse Mercator (South Orientated)',27), ...
         ValuePair('Unknown',32767)};

