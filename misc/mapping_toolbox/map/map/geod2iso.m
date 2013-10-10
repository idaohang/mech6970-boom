function latout = geod2iso(varargin)
%GEOD2ISO  Convert geodetic latitude to isometric latitude.
%   GEOD2ISO is obsolete; use CONVERTLAT.
%
%   lat = GEOD2ISO(lat0) computes the isometric latitude given the
%   geodetic latitude, using the default Earth ellipsoid from ALMANAC.
%
%   lat = GEOD2ISO(lat0,ELLIPSOID) uses the ellipsoid definition given in
%   the input vector ellipsoid.  ELLIPSOID can be determined from the ALMANAC
%   function.  If omitted, the default Earth ellipsoid is assumed.
%
%   lat = GEOD2ISO(lat0,'units') uses the units defined by the input string
%   'units'.  If omitted, default units of degrees are assumed.
%
%   lat = GEOD2ISO(lat0,ELLIPSOID,'units') uses the ellipsoid and 'units'
%   definitions provided by the corresponding inputs.

%  Copyright 1996-2007 The MathWorks, Inc.
%  $Revision: 1.11.4.2 $  $Date: 2007/02/11 05:47:01 $

latout = doLatitudeConversion(mfilename,'geodetic','isometric',varargin{:});
