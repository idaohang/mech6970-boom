function latout = aut2geod(varargin)
%AUT2GEOD  Convert authalic latitude to geodetic latitude.
%   AUT2GEOD is obsolete; use CONVERTLAT.
%
%   lat = AUT2GEOD(lat0) converts from the authalic latitude to the
%   geodetic latitude, using the default Earth ellipsoid from ALMANAC.
%
%   lat = AUT2GEOD(lat0,ELLIPSOID) uses the ellipsoid definition given in
%   the input vector ellipsoid.  ELLIPSOID can be determined from the
%   ALMANAC function.
%
%   lat = AUT2GEOD(lat0,'units') uses the units defined by the input string
%   'units'.  If omitted, default units of degrees are assumed.
%
%   lat = AUT2GEOD(lat0,ELLIPSOID,'units') uses the ellipsoid and 'units'
%   definitions provided by the corresponding inputs.

%  Copyright 1996-2007 The MathWorks, Inc.
%  $Revision: 1.10.4.2 $  $Date: 2007/02/11 05:46:48 $

latout = doLatitudeConversion(mfilename,'authalic','geodetic',varargin{:});
