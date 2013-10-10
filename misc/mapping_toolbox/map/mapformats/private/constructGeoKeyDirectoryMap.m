function [idToNameMap, nameToIdMap] = constructGeoKeyDirectoryMap
%CONSTRUCTGEOKEYDIRECTORYMAP Construct maps for GeoKeyDirectoryTag
%
%   [idToNameMap, nameToIdMap] = constructGeoKeyDirectoryMap maps GeoTIFF
%   key ID (GeoKey) to key name. The function returns one or two maps based
%   on the number of requested outputs. If one output is requested, the
%   output parameter maps the GeoKey ID to its name. If two outputs are
%   requested, the second output parameter maps the key name (in lower
%   case) back to the GeoKey ID. The names and values are obtained from the
%   GeoTIFF specification.

% Copyright 2010 The MathWorks, Inc.
% $Revision: 1.1.6.1 $  $Date: 2010/09/24 14:33:58 $

keyName = { ...
    'GTModelTypeGeoKey', ...
    'GTRasterTypeGeoKey', ...
    'GTCitationGeoKey', ...
    'GeographicTypeGeoKey', ...
    'GeogCitationGeoKey', ...           
    'GeogGeodeticDatumGeoKey', ...      
    'GeogPrimeMeridianGeoKey', ...      
    'GeogLinearUnitsGeoKey', ...        
    'GeogLinearUnitSizeGeoKey', ...     
    'GeogAngularUnitsGeoKey', ...       
    'GeogAngularUnitSizeGeoKey', ...    
    'GeogEllipsoidGeoKey', ...          
    'GeogSemiMajorAxisGeoKey', ...      
    'GeogSemiMinorAxisGeoKey', ...     
    'GeogInvFlatteningGeoKey', ...      
    'GeogAzimuthUnitsGeoKey', ...       
    'GeogPrimeMeridianLongGeoKey', ...  
    'ProjectedCSTypeGeoKey', ...          
    'PCSCitationGeoKey', ...              
    'ProjectionGeoKey', ...               
    'ProjCoordTransGeoKey', ...           
    'ProjLinearUnitsGeoKey', ...          
    'ProjLinearUnitSizeGeoKey', ...       
    'ProjStdParallel1GeoKey', ...         
    'ProjStdParallel2GeoKey', ...         
    'ProjNatOriginLongGeoKey', ...       
    'ProjNatOriginLatGeoKey', ...         
    'ProjFalseEastingGeoKey', ...         
    'ProjFalseNorthingGeoKey', ...        
    'ProjFalseOriginLongGeoKey', ...
    'ProjFalseOriginLatGeoKey', ...
    'ProjFalseOriginEastingGeoKey', ...   
    'ProjFalseOriginNorthingGeoKey', ...  
    'ProjCenterLongGeoKey', ...           
    'ProjCenterLatGeoKey', ...            
    'ProjCenterEastingGeoKey', ...        
    'ProjCenterNorthingGeoKey', ...       
    'ProjScaleAtNatOriginGeoKey', ...     
    'ProjScaleAtCenterGeoKey', ...        
    'ProjAzimuthAngleGeoKey', ...         
    'ProjStraightVertPoleLongGeoKey', ... 
    'VerticalCSTypeGeoKey', ...           
    'VerticalCitationGeoKey', ...         
    'VerticalDatumGeoKey', ...            
    'VerticalUnitsGeoKey', ...
    'UserDefined'};

keyID = { ...
   1024, ... 
   1025, ...
   1026, ...
   2048, ... 
   2049, ... 
   2050, ... 
   2051, ... 
   2052, ... 
   2053, ...
   2054, ...
   2055, ... 
   2056, ... 
   2057, ... 
   2058, ... 
   2059, ... 
   2060, ... 
   2061, ... 
   3072, ... 
   3073, ... 
   3074, ... 
   3075, ... 
   3076, ...
   3077, ... 
   3078, ... 
   3079, ... 
   3080, ... 
   3081, ... 
   3082, ... 
   3083, ... 
   3084, ... 
   3085, ... 
   3086, ... 
   3087, ...
   3088, ... 
   3089, ... 
   3090, ...
   3091, ... 
   3092, ...
   3093, ... 
   3094, ... 
   3095, ... 
   4096, ... 
   4097, ... 
   4098, ... 
   4099, ...
   32767};

idToNameMap = containers.Map(keyID, keyName);
if nargout > 1
    nameToIdMap = containers.Map(lower(keyName), keyID);
end
