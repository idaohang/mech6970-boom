function corner = firstCorner(center, delta)
% Working in a single dimension of a rectilinear, spatially-referenced
% image or raster data set, given the coordinate of the center of the 1,1
% cell, CENTER, and the signed cell width or height, DELTA, compute the
% coordinate of the outer corner in the same dimension. Perform the
% computation in a robust manner, working with whole numbers a much as
% possible.

% Copyright 2010 The MathWorks, Inc.
% $Revision: 1.1.6.1 $  $Date: 2010/10/11 14:47:40 $

wholeNumberOffset = fix(center);

[oNum, oDen] = rat(center - wholeNumberOffset);
[dNum, dDen] = rat(delta);

if abs(wholeNumberOffset + (oNum / oDen) - center) < 5*abs(center) ...
        && delta == dNum / dDen
    % Rational approximation is reasonable with respect to machine precision
    
    % Least common multiple of center denominator and 2x delta denominator
    lcmDen = lcm(oDen, 2 * dDen);
    
    corner = wholeNumberOffset ...
        + (2 * oNum * (lcmDen/oDen) + dNum * (lcmDen/dDen)) / (2 * lcmDen);
else
    % Rational approximation is worse than machine precision
    corner = center + delta - delta / 2;
end
