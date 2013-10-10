function rotatetext(h,direction)
%ROTATETEXT Rotate text to projected graticule
% 
%   ROTATETEXT rotates displayed text objects to account for the
%   curvature of the graticule.  The objects are selected interactively
%   from a graphical user interface.
%   
%   ROTATETEXT(OBJECTS) rotates the selected objects.  OBJECTS may be a
%   name string recognized by HANDLEM, or a vector of handles to
%   displayed text objects.
%   
%   ROTATETEXT(OBJECTS,DIRECTION) accepts a DIRECTION string, which must
%   be either 'forward' (the default value) or 'inverse'. Specifying
%   'inverse' causes ROTATETEXT to remove the rotation added by an
%   earlier call to ROTATETEXT.
%
%   Meridian and parallel labels can be rotated automatically by setting
%   the map axes LabelRotation property to on. ROTATETEXT does not
%   support the 'globe' projection.

% Copyright 1996-2009 The MathWorks, Inc.
% $Revision: 1.5.4.5 $  $Date: 2009/12/02 06:43:43 $

error(nargchk(0, 2, nargin, 'struct'))

% Get handles to displayed text objects
if nargin < 1
    h = handlem;
elseif ischar(h)
    h = handlem(h);
elseif ~ishghandle(h);
	error('map:rotatetext:invalidObject', ...
        'Object must be a name string or handle.')
end

% Get direction / set forward to true or false
if nargin < 2
    forward = true;
else
    switch(direction)
        case 'forward'
            forward = true;
        case 'inverse'
            forward = false;
        otherwise
            error('map:rotatetext:invalidDirectionString', ...
                '%s must be ''%s'' or ''%s''.', 'DIRECTION','forward','inverse')
    end
end

% Validate map axes
mstruct = gcm;
if strcmp(mstruct.mapprojection,'globe')
    error('map:rotatetext:usingGlobe', ...
        '%s does not work with the globe projection.','ROTATETEXT')
end

% open limits to avoid bumping against the frame

t = defaultm(mstruct.mapprojection); % temporary mstruct, in degrees
mstruct.flatlimit = fromDegrees(mstruct.angleunits, t.trimlat);
mstruct.flonlimit = fromDegrees(mstruct.angleunits, t.trimlon);

% calculate vector rotation introduced by the projection

for i=1:length(h)
    if ishghandle(h(i),'text')
        pos = get(h(i),'position');
        [lat,lon] = minvtran(pos(1),pos(2));
        th = vfwdtran(mstruct,lat,lon,fromDegrees(mstruct.angleunits,90));
        rot = get(h(i),'Rotation');
        if forward
            set(h(i),'rotation',th+rot);
        else
            set(h(i),'rotation',-th+rot);
        end
    end
end

%-----------------------------------------------------------------------

function th = vfwdtran(mstruct, lat, lon, azf)
%VFWDTRAN  Direction angle in map plane from azimuth on sphere
%
%   th = VFWDTRAN(MSTRUCT, LAT, LON, AZF) transforms the azimuth angle
%   AZF at the specified point on the sphere into the projected map
%   coordinate system defined by MSTRUCT. Inputs LAT and LON are scalar
%   latitude and longitude. LAT, LON, and AZF are in the angle units
%   specified by MSTRUCT. The output angle TH is measured
%   counterclockwise from the X-axis, and is always in degrees. If the
%   input (LAT, LON) falls outside the map frame or if we are unable to
%   compute a vector connecting it to another nearby point, then a value
%   of 0 is returned.
%
%   This is a specialized version of the public VFWDTRAN function.

mapproj = mstruct.mapprojection;
[x, y] = feval(mapproj, mstruct, lat,  lon,  'geopoint', 'forward');
if isempty(x)
    % Point falls outside the map frame; leave the text unrotated.
    th = 0;
    return
end

% Compute points about 10 centimeters away in both directions. Go in
% both directions in case the point (LAT, LON) is on or close to the
% frame boundary. Assume a sphere when using reckon; it's simpler and
% provides more than enough precision for text rotation.
angleUnits = mstruct.angleunits;
rng = 10*epsm(angleUnits);
azb = azf + fromDegrees(angleUnits,180);
[latf, lonf] = reckon('rh', lat, lon, rng, azf, angleUnits);
[latb, lonb] = reckon('rh', lat, lon, rng, azb, angleUnits);

% Project the nearby points. Note that if (lat,lon) falls near the edge
% of the map frame, then these may be trimmed away or they may project
% to a location on the far side of the map.
[xf, yf] = feval(mapproj, mstruct, latf, lonf, 'geopoint', 'forward');
[xb, yb] = feval(mapproj, mstruct, latb, lonb, 'geopoint', 'forward');

% Vector lengths in map plane.
if ~isempty(xf)
    df = hypot(yf - y, xf - x);
else
    df = Inf;
end

if ~isempty(xb)
    db = hypot(yb - y, xb - x);
else
    db = Inf;
end

% Convert vectors in the plane to angles in degrees; if there's a
% choice, use the result from the shortest vector.
th = 0;
if isfinite(df)
    th = radtodeg(atan2(yf - y, xf - x));
end
if db < df
    th = radtodeg(atan2(y - yb, x - xb));
end
