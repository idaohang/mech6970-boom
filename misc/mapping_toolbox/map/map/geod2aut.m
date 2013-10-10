function latout = geod2aut(varargin)
%GEOD2AUT  Convert geodetic latitude to authalic latitude.
%   GEOD2AUT is obsolete; use CONVERTLAT.
%
%   lat = GEOD2AUT(lat0) converts from the geodetic latitude to the
%   authalic latitude, using the default Earth ellipsoid from ALMANAC.
%
%   lat = GEOD2AUT(lat0,ELLIPSOID) uses the ellipsoid definition given in
%   the input vector ellipsoid.  ELLIPSOID can be determined from the
%   ALMANAC function.
%
%   lat = GEOD2AUT(lat0,'units') uses the units defined by the input string
%   'units'.  If omitted, default units of degrees are assumed.
%
%   lat = GEOD2AUT(lat0,ELLIPSOID,'units') uses the ellipsoid and 'units'
%   definitions provided by the corresponding inputs.

%  Copyright 1996-2007 The MathWorks, Inc.
%  $Revision: 1.10.4.2 $  $Date: 2007/02/11 05:46:58 $

latout = doLatitudeConversion(mfilename,'geodetic','authalic',varargin{:});
