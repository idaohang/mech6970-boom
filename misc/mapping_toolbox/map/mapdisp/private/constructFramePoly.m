function [latfrm, lonfrm] = constructFramePoly(mstruct)

% Copyright 2009 The MathWorks, Inc.
% $Revision: 1.1.6.3 $  $Date: 2009/03/30 23:39:23 $

framelat  = mstruct.flatlimit;
framelon  = mstruct.flonlimit;
units     = mstruct.angleunits;
fillpts   = mstruct.ffill;

epsilon   = 1000*epsm('degrees');

if all(~isinf(framelat))
    % Frame is not an azimuthal disk

    framelat = sort(toDegrees(units,framelat));  % Convert input
    framelon = sort(toDegrees(units,framelon));  % data to degrees

    framelat = framelat + [epsilon -epsilon];  %  Avoid clipping at edge of
    framelon = framelon + [epsilon -epsilon];  %  of map

    lats = linspace(min(framelat),max(framelat),fillpts)';   % Fill vectors with
    lons = linspace(min(framelon),max(framelon),fillpts)';   % frame limits

    latfrm = [lats;           framelat(2)*ones(size(lats));  % Construct
              flipud(lats);   framelat(1)*ones(size(lats))]; % complete frame
    lonfrm = [framelon(1)*ones(size(lons));    lons;         % vectors
              framelon(2)*ones(size(lons));    flipud(lons);];
          
elseif any(isinf(framelat))
    %  Frame is an azimuthal disk

    if strcmp(mstruct.mapprojection,'vperspec')
        % vertical perspective requires special treatment
        P = mstruct.mapparallels/mstruct.geoid(1) + 1;
        % reset the frame
        trimlat = [-inf  min( [ acos(1/P)-5*epsm('radians')  max(framelat)  1.5533] ) ]; % 1.5533 rad = 89 degrees
        mstruct.flatlimit = radtodeg(trimlat);
        framelat  = mstruct.flatlimit;
    end

    framelat = toDegrees(units,max(framelat)); % Convert disk radius
    framelat = framelat - epsilon;  % Avoid clipping at map edge

    az = [0 360];      % Compute azimuthal points on frame

    [latfrm,lonfrm] = scircle1('gc',0,0,framelat,az,[],'degrees',fillpts);
end

% Work in degrees, but return output in the angleunits of the map axes
[latfrm,lonfrm] = fromDegrees(mstruct.angleunits,latfrm,lonfrm);
