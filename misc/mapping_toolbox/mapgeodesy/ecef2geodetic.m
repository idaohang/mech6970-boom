function [phi, lambda, h] = ecef2geodetic(x, y, z, ellipsoid)
%ECEF2GEODETIC Convert geocentric (ECEF) to geodetic coordinates
%
%   [PHI, LAMBDA, H] = ECEF2GEODETIC(X, Y, Z, ELLIPSOID) converts point
%   locations in geocentric Cartesian coordinates, stored in the
%   coordinate arrays X, Y, Z, to geodetic coordinates PHI (geodetic
%   latitude in radians), LAMBDA (longitude in radians), and H (height
%   above the ellipsoid). The geodetic coordinates refer to the
%   reference ellipsoid specified by ELLIPSOID (a row vector with the
%   form [semimajor axis, eccentricity]). X, Y, and Z must use the same
%   units as the semimajor axis;  H will also be expressed in these
%   units.  X, Y, and Z must have the same shape; PHI, LAMBDA, and H
%   will have this shape also.
%
%   For a definition of the geocentric system, also known as
%   Earth-Centered, Earth-Fixed (ECEF), see the help for GEODETIC2ECEF.
%
%   See also ECEF2LV, GEODETIC2ECEF, GEOCENTRIC2GEODETICLAT, LV2ECEF.

% Copyright 2005-2009 The MathWorks, Inc.
% $Revision: 1.1.6.4 $  $Date: 2009/04/15 23:34:43 $

% Reference
% ---------
% Paul R. Wolf and Bon A. Dewitt, "Elements of Photogrammetry with
% Applications in GIS," 3rd Ed., McGraw-Hill, 2000 (Appendix F-3).

% Implementation Notes from Rob Comer
% -----------------------------------
% The implementation below follows Wolf and DeWitt quite literally,
% with a few important exceptions required to ensure good numerical
% behavior:
%
% 1) I used ATAN2 rather than ATAN in the formulas for beta and phi.  This
%    avoids division by zero (or a very small number) for points on (or
%    near) the Z-axis.
%
% 2) Likewise, I used ATAN2 instead of ATAN when computing beta from phi
%    (conversion from geodetic to parametric latitude), ensuring
%    stability even for points at very high latitudes.
%
% 3) Finally, I avoided dividing by cos(phi) -- also problematic at high
%    latitudes -- in the calculation of h, the height above the ellipsoid.
%    Wold and Dewitt give
%
%                   h = sqrt(X^2 + Y^2)/cos(phi) - N.
%
%    The trick is to notice an alternative formula that involves division
%    by sin(phi) instead of cos(phi), then take a linear combination of the
%    two formulas weighted by cos(phi)^2 and sin(phi)^2, respectively. This
%    eliminates all divisions and, because of the identity cos(phi)^2 +
%    sin(phi)^2 = 1 and the fact that both formulas give the same h, the
%    linear combination is also equal to h.
%
%    To obtain the alternative formula, we simply rearrange
%
%                   Z = [N(1 - e^2) + h]sin(phi)
%    into
%                   h = Z/sin(phi) - N(1 - e^2).
%
%    The linear combination is thus
%
%        h = (sqrt(X^2 + Y^2)/cos(phi) - N) cos^2(phi)
%            + (Z/sin(phi) - N(1 - e^2))sin^2(phi)
%
%    which simplifies to
%
%      h = sqrt(X^2 + Y^2)cos(phi) + Zsin(phi) - N(1 - e^2sin^2(phi)).
%
%    From here it's not hard to verify that along the Z-axis we have
%    h = Z - b and in the equatorial plane we have h = sqrt(X^2 + Y^2) - a.

% Ellipsoid constants
a  = ellipsoid(1);       % Semimajor axis
e2 = ellipsoid(2) ^ 2;   % Square of first eccentricity
ep2 = e2 / (1 - e2);     % Square of second eccentricity
f = 1 - sqrt(1 - e2);    % Flattening
b = a * (1 - f);         % Semiminor axis

% Longitude
lambda = atan2(y,x);

% Distance from Z-axis
rho = hypot(x,y);

% Bowring's formula for initial parametric (beta) and geodetic (phi) latitudes
beta = atan2(z, (1 - f) * rho);
phi = atan2(z   + b * ep2 * sin(beta).^3,...
            rho - a * e2  * cos(beta).^3);

% Fixed-point iteration with Bowring's formula
% (typically converges within two or three iterations)
betaNew = atan2((1 - f)*sin(phi), cos(phi));
count = 0;
while any(beta(:) ~= betaNew(:)) && count < 5
    beta = betaNew;
    phi = atan2(z   + b * ep2 * sin(beta).^3,...
                rho - a * e2  * cos(beta).^3);
    betaNew = atan2((1 - f)*sin(phi), cos(phi));
    count = count + 1;
end

% Calculate ellipsoidal height from the final value for latitude
sinphi = sin(phi);
N = a ./ sqrt(1 - e2 * sinphi.^2);
h = rho .* cos(phi) + (z + e2 * N .* sinphi) .* sinphi - N;
