function [latTrimmed, lonTrimmed] = ...
    trimPolygonToSmallCircle(lat, lon, latCenter, lonCenter, radius, inc)
%TRIMPOLYGONTOSMALLCIRCLE Trim lat-lon polygon to small circle
%
%   [latTrimmed, lonTrimmed] = ...
%       trimPolygonToSmallCircle(lat, lon, latCenter, lonCenter, radius, inc)
%   trims the polygon defined by vectors LAT and LON to the small circle
%   defined by (latCenter, lonCenter) and RADIUS.  LAT and LON may contain
%   multiple rings separated by NaNs.  Outer rings should be clockwise,
%   inner rings counterclockwise.  INC is the distance increment to be used
%   to define the edges of trimmed polygons that intersect edge of the
%   small circle itself.  All inputs and outputs are assumed to be in units
%   of radians.

% Copyright 2005-2010 The MathWorks, Inc.
% $Revision: 1.1.8.9 $  $Date: 2010/06/07 16:34:24 $

% Keep track of the shape of the input vectors, because polyjoin will
% automatically convert everything to column vectors.
rowVectorInput = (size(lat,2) > 1);
lat = lat(:);
lon = lon(:);

% Make sure lat and lon arrays are NaN-terminated.
nanTerminatedInput = isnan(lon(end));
if ~nanTerminatedInput
    lon(end+1,1) = NaN;
    lat(end+1,1) = NaN;
end

% Tolerance for snapping to limits.
tolSnap = 10*eps(pi);

% Set the tolerance for closing nearly-closed rings to 5 degrees.
tolClose = 5 * pi/180;

% Pre-process lat-lon polygons.
[lat, lon] = preprocessLatLonPolygons(lat, lon, tolSnap, tolClose);

% Make sure that radius is slightly less than pi.  This buffer allows us
% to construct a ring around the antipode, if necessary. Set the buffer
% width to 0.5 degrees.
bufferWidth = 0.5 * (pi/180);
if radius > (pi - bufferWidth)
    radius = radius - bufferWidth;
end

% Transform latitude-longitude to range-azimuth.
[rng, az] = latlon2rngaz(lat, lon, latCenter, lonCenter);

% Save a copy of the input polygon.
rngIn = rng;
azIn = az;

% Trim range-azimuth vectors such that rng <= radius, leaving "dangling
% ends" along the circumference of the circle.
[rng, az] = truncateAtBoundary(rng, az, radius, false);
if all(isnan(rng))
    rng = [];
    az = [];
end

% When trimming a topologically-consistent, multi-part polygon to a
% small circle centered on the projection origin (the "trimming
% circle"), we need to determine when, in order to properly enclose all
% interior areas in the projected map plane, it's necessary to add a
% ring that traces the trimming circle and encloses the origin. There
% are four such cases. Two are handled by subfunction originFallsInside
% and the other two are handled by trimmingCircleFallsInside.
if ~isempty(rng)
    % Trace curves that have been disconnected by the truncation process,
    % interpolating additional vertices along the circumference.
    [rng, az, closureNeeded] = closePolygonInCircle(rng, az, radius, inc);

    % Check to see if there's a ring (after the trimming step above)
    % that encloses the trimming circle itself. This check is necessary
    % only if function closePolygonInCircle has not already closed a
    % ring that touches the trimming circle.
    if ~closureNeeded && trimmingCircleFallsInside(rng, az)
        % All the "outermost" rings are counter-clockwise and thus
        % enclose the trimming circle itself; enclose them within a
        % clockwise ring tracing the circumference of the circle.
        [rngCircle, azCircle] = smallCircle(radius, inc);
        rng = [rng; NaN; rngCircle];
        az  = [ az; NaN; azCircle];
    end
else
    % All the inputs have been trimmed away, but check the full set of data
    % to see if the origin (along with the entire interior of the trimming
    % circle) falls on the inside of the input polygon.
    if originFallsInside(rngIn, azIn)
        % Return a clockwise ring tracing the circumference of the circle.
        [rng, az] = smallCircle(radius, inc);
    end
end

% Make sure terminating NaNs haven't been lost.
if ~isempty(rng) && ~isnan(rng(end))
    rng(end+1,1) = NaN;
    az(end+1,1) = NaN;
end

% Transform range-azimuth back to latitude-longitude.
[latTrimmed, lonTrimmed] = rngaz2latlon(rng, az, latCenter, lonCenter);

% Restore shape if necessary.
if rowVectorInput
    latTrimmed = latTrimmed';
    lonTrimmed = lonTrimmed';
end

%-----------------------------------------------------------------------

function [rng, az] = smallCircle(radius, inc)
% Return a ring enclosing a small circle with the specified radius.

n = ceil(2*pi/inc);
az = 2*pi*(0:n)'/n;
rng = radius + zeros(size(az));

%--------------------------------------------------------------------------

function [rng, az, closureNeeded] = closePolygonInCircle(rng, az, radius, inc)
% Trace and re-connect open curves which start and/or end at the
% trimming radius. Be sure to work in azimuth modulo 2*pi.

[first,last] = internal.map.findFirstLastNonNan(az);

% Identify open curves that start and/or end at rng == radius, allowing
% for round off of azimuth values during a previous unwrapping step.
isOpen = (rng(first) == radius | rng(last) == radius) ...
    & ~(rng(first) == rng(last)  ...
        & abs(wrapToPi(az(first) - az(last))) < 100*eps(pi));

closureNeeded = any(isOpen);
if closureNeeded
    % Number of vertices needed to interpolate a full circle
    nCircumference = ceil(2*pi/inc);
    
    % Allocate output arrays, allowing for additional vertices.
    rngTraced = NaN(numel(rng) + nCircumference,1);
    azTraced = rngTraced;
    
    % Indices for open curves
    firstOpen = first(isOpen);
    lastOpen  = last(isOpen);

    % Construct a lookup table which, given an open curve index k
    % returns the index of the open curve next(k) whose start point
    % coincides with or follows (in terms of increasing azimuth, modulo
    % 2*pi) curve k's end point.
    next = nextCurveLookup(rng, az, radius, firstOpen, lastOpen);
    
    % Trace the open curves, copying vertices into rngTraced and
    % azTraced.  n, a positive, scalar integer indicates where to start
    % when copying addition vertices (from closed curves) into rngTraced
    % and azTraced.
    [~,azCircle] = smallCircle(radius, inc);
    [rngTraced, azTraced, n] = traceOpenCurves(rng, az, ...
        firstOpen, lastOpen, next, rngTraced, azTraced, radius, azCircle);
    
    % Append curves that were closed already.
    [rngTraced, azTraced] = appendClosedCurves(rng, az, ...
        first(~isOpen), last(~isOpen), rngTraced, azTraced, n + 1);  
    
    % Correct for any excess allocation.
    [rng, az] = removeExtraNanSeparators(rngTraced, azTraced);
    
    % Wrap the azimuths.
    az = mod(az, 2*pi);

    % Remove duplicate vertices.
    duplicate = [diff(az)==0 & diff(rng)==0; false];
    az(duplicate) = [];
    rng(duplicate) = [];
end

%--------------------------------------------------------------------------

function next = nextCurveLookup(rng, az, radius, firstOpen, lastOpen)
% Construct a lookup table which, given an open curve index k
% returns the index of the open curve next(k) whose start point
% coincides with or follows (in terms of increasing azimuth, modulo
% 2*pi) curve k's end point.

tolAzimuth = 100*eps(pi);
tolRange = eps(10);

% Start points as column vectors
rs = rng(firstOpen);
as = az(firstOpen);

% End points as row vectors
re = rng(lastOpen)';
ae = az(lastOpen)';

endsOnCircle   = (re == radius);
startsOnCircle = (rs == radius);

aeCircle = ae(endsOnCircle);
asCircle = as(startsOnCircle);

% Expand into n-by-n arrays.
n = numel(firstOpen);
nOnes = ones(1,n);

rs = rs(:,nOnes);
as = as(:,nOnes);

re = re(nOnes,:);
ae = ae(nOnes,:);

% n-by-n logical connects is true when an end point coincides
% (or nearly coincides) with a start point.
connects = (abs(wrapToPi(as - ae)) < tolAzimuth) ...
    & ((re == rs) | (re < radius & rs < radius & abs(re-rs) < tolRange));

% Account also for the case in which the start and end points are
% distinct, but a start point on the circle immediately follows an end
% point when the circle is traversed in a clockwise direction.
n = numel(aeCircle);
nOnes = ones(1,n);
asCircle = asCircle(:,nOnes);
aeCircle = aeCircle(nOnes,:);
[~,nextOnCircle] = min(mod(asCircle - aeCircle, 2*pi));
connectsOnCircle = false(n,n);
for k = 1:numel(nextOnCircle)
    connectsOnCircle(nextOnCircle(k),k) = true;
end

% Combine results for both types of connectivity.
connects(startsOnCircle, endsOnCircle) = connectsOnCircle;
[next, ~] = find(connects);

%--------------------------------------------------------------------------

function [rLinked, aLinked, n] = traceOpenCurves(rng, az, ...
    firstOpen, lastOpen, next, rLinked, aLinked, radius, azCircle)
% Trace to link up open curves and add vertices on circumference of circle.

traced = false(size(firstOpen));
nOpen = numel(traced);
nTraced = 0;
k = 1; % Index to current open curve
n = 1; % Index to current position in output vertex arrays
f = 1; % Index to start of current curve in set of linked curves
while any(~traced)
    nTraced = nTraced + 1;
    assert(nTraced <= nOpen, ...
        'map:trimPolygonToSmallCircle:tracingFailed', ...
        'Failed to converge when tracing open curves.')
            
    if traced(k)
        k = find(~traced);
        k = k(1);
        n = n + 1; % Allow for NaN-separator
    end
    s = firstOpen(k);
    e = lastOpen(k);
    m = n + e - s;
    rLinked(n:m) = rng(s:e);
    aLinked(n:m) = az(s:e);
    traced(k) = true;
    k1 = k;
    k = next(k);
    if (rng(lastOpen(k1)) == radius) && (rng(firstOpen(k)) == radius)
        % Curves end and start on circle:
        %   Insert extra vertices along the circumference.
        az1 = mod(az(lastOpen(k1)),2*pi);
        az2 = mod(az(firstOpen(k)),2*pi);
        if az1 <= az2
            indx = find(az1 < azCircle & azCircle < az2);
        else
            indx = [find(az1 < azCircle); find(azCircle < az2)];
        end
        n = m + 1;
        m = n + numel(indx) - 1;
        rLinked(n:m) = radius;
        aLinked(n:m) = azCircle(indx);
    end
    if traced(k)
        % Close up curve and replicate first vertex if needed.
        m = m + 1;
        rLinked(m) = rLinked(f);
        aLinked(m) = aLinked(f);
        f = m + 2;
    end
    n = m + 1;
end

%--------------------------------------------------------------------------

function [u, v, n] = appendClosedCurves( ...
    x, y, firstClosed, lastClosed, u, v, n)
% Copy closed curves from NaN-separated vertex arrays (x,y) to
% NaN-separated vertex arrays (u,v). The positive integer n is the index
% of a starting vertex in (u,v), n.

% Adapted from toolbox/map/map/private/gluePolygonsOnVerticalEdges.m

for k = 1:numel(firstClosed)
    % First and last indices of k-th closed curve in (x,y).
    s = firstClosed(k);
    e = lastClosed(k);
    
    % Compute index m for (u,v) such that n:m is the same size as s:e
    m = n + e - s;
    
    % Copy vertices from k-th closed curve.
    u(n:m) = x(s:e);
    v(n:m) = y(s:e);
    
    % Advance by 2 instead of 1, leaving a NaN-separator in u and v to
    % separate this curve from the next one.
    n = m + 2;
end

%--------------------------------------------------------------------------

function tf = originFallsInside(rng, az)
% Analyze the situation in which no polygon vertices at all fall within
% the trimming circle and either:
% 
% 1a. There is a ring that wraps the origin in a planar azimuthal system
% centered on the origin and encloses the origin (the origin falls on
% its right-hand side; equivalently, the ring is CW in the planar system;
% also, the ring encloses not just the origin but rather all points
% within the trimming circle), and this ring does not enclose any other
% ring that wraps the origin but does not enclose it (is CCW in the
% planar system and excludes all points within the trimming circle), or
% 
% 1b. No rings at all wrap and enclose the origin, but there is a ring
% (CCW in the planar system) that encloses the origin (or any point
% within the trimming circle) because it is not itself enclosed by
% another ring (CW in the planar system) that does not enclose the
% origin (and thus excludes the origin and all points within the
% trimming circle).

az = snapNearlyEqualAzimuths(az);
[cw, wrapsOrigin] = ispolycwPolar(-az,rng);
if any(wrapsOrigin & cw)
    % At least one CW ring wraps the origin. Output true if the ring
    % that most tightly wraps the origin (the origin-wrapping ring with
    % the smallest area, that is) is CW (case 1a).    
    tf = cw(findPoleWrappingRing(rng, az, wrapsOrigin, @min));
elseif any(~wrapsOrigin & ~cw)
    % There are one or more CCW rings that do not wrap the origin.
    % Output true if at least one CCW ring is not enclosed by a CW ring
    % -- which may or may not wrap the origin (case 1b).
    tf = ccwRingsAreUnenclosed(rng, az);
else
    % If there are any CCW rings, all of them wrap the origin.
    tf = false;
end

%--------------------------------------------------------------------------

function tf = trimmingCircleFallsInside(rng, az)
% Analyze the situation in which at least one ring falls within the
% trimming circle, but no ring actually intersects the trimming circle,
% and, considering only the rings defined by vertices that fall within
% the trimming circle, there is at least one ring that encloses the
% trimming circle itself and is unenclosed by another other ring. This
% can happen in one of two ways, either:
%
% 2a. There is a ring that wraps the origin in a planar azimuthal system
% centered on the origin and excludes the origin (the origin falls on
% its left-hand side; equivalently, the ring is CCW in the planar
% system) and this ring is not enclosed by another other ring (which, if
% it existed, would have to wrap the origin and be CW in the planar
% system), or
%
% 2b. There is at least one ring that does not wrap the origin but
% encloses the trimming circle (and which must be CCW in the planar
% system) and which it is not itself enclosed by a ring that does not
% also enclose the trimming circle. (Such a ring must be CW in the
% planar system.)

az = snapNearlyEqualAzimuths(az);
[cw, wrapsOrigin] = ispolycwPolar(-az,rng);
if any(wrapsOrigin & ~cw)
    % At least one CCW ring wraps the origin. Output true if the ring
    % that most loosely wraps the origin (the origin-wrapping ring with
    % the largest area, that is) is CCW (case 2a).
     tf = ~cw(findPoleWrappingRing(rng, az, wrapsOrigin, @max));
else
    % There may be one or more CCW rings, but none of them wrap the
    % origin. Output true if there is at least one CCW ring that is not
    % enclosed by a CW ring -- which may or may not wrap the origin
    % (case 2b).
    tf = ccwRingsAreUnenclosed(rng, az);
end

%--------------------------------------------------------------------------

function k = findPoleWrappingRing(rng, az, wrapsOrigin, extremum)
% Find the tightest pole-wrapping ring if extremum is @min or find the
% loosest pole-wrapping ring if extremum is @max.

% Locate the rings within the rng and az vectors.
[first, last] = internal.map.findFirstLastNonNan(rng);

% Compute the areas of the origin-wrapping rings, and ignore rings
% that don't wrap the origin by setting their areas to Inf.
if isequal(extremum,@min)
    ringArea = Inf(size(first));
else
    ringArea = zeros(size(first));
end

for k = find(wrapsOrigin(:))' % Make sure right hand side is row vector
    ringArea(k) = abs(areaPolar(...
        az(first(k):last(k)), rng(first(k):last(k))));
end

% Find the index of the minimum or maximum-area origin-wrapping ring. We
% expect the result to always be scalar because any(wrapsOrigin) will
% always be true and rings must have unique areas because they are not
% allowed to intersect.
k = find(ringArea == extremum(ringArea));

%--------------------------------------------------------------------------

function tf = ccwRingsAreUnenclosed(rng, az)
% Working in range-azimuth space, return false if and only if every
% counter-clockwise ring is enclosed within a clockwise ring. In other
% words, return true if any counter-clockwise ring is unenclosed.
%
% Note that one _cannot_ simply apply pol2cart and then work in a
% Cartesian system. The reason is that a linear segment between two
% vertices that is a straight line in range-azimuth will map to a curve
% in the Cartesian system, so to do that we'd have to first interpolate
% extra vertices and select a threshold for a sufficiently dense
% sampling, etc., and all results would still be approximate.

az = snapNearlyEqualAzimuths(az);
cw = ispolycwPolar(-az,rng);
if any(~cw)
    if all(~cw)
        % All rings are counter-clockwise (CCW).
        tf = true;
    else
        % There's a combination of clockwise and counter-clockwise rings.
        
        % Locate the rings within the rng and az vectors.
        [first, last] = internal.map.findFirstLastNonNan(rng);
        
        % For test points use the first element of each CCW ring.
        rngTest = rng(first(~cw));
        azTest  = az( first(~cw));
        
        % Initialize a logical vector to track the status of each CCW ring.
        enclosed = false(size(azTest));
        
        % Isolate the clockwise rings.
        firstCW = first(cw);
        lastCW = last(cw);
        
        for k = 1:numel(firstCW)
            % See which test points are enclosed by the k-th CW polygon
            % and update status according.
            rngCW = rng(firstCW(k):lastCW(k));
            azCW  = az( firstCW(k):lastCW(k));
            enclosed = enclosed | inSingleRing(azTest, rngTest, azCW, rngCW);
        end
        
        tf = any(~enclosed);
    end
else
    % There are no counter-clockwise rings at all.
    tf = false;
end

%-----------------------------------------------------------------------

function az = snapNearlyEqualAzimuths(az)

[first, last] = internal.map.findFirstLastNonNan(az);
q = (abs(az(last) - az(first)) <= eps(4*pi));
az(last(q)) = az(first(q));

%-----------------------------------------------------------------------

function tf = inSingleRing(azTest, rngTest, az, rng)
% Return true if and only if the points defined by azTest and rngTest
% fall within (or on the edge of) the ring defined by the vertex arrays az
% and rng. az and rng define a single ring, without NaN-separators or
% terminating NaNs. The ring may or may not wrap the origin.

% The basic strategy is to convert the ring to a simple, closed polygon
% in azimuth-range coordinates by trimming it to an interval of width
% 2*pi in azimuth.

wrapsOrigin = abs(az(end) - az(1)) > pi;
if wrapsOrigin
    % If the ring wraps the origin, so that az(end) is equal (or nearly
    % equal) to az(1), keep the computations somewhat simpler by
    % choosing the interval to match the start point of the ring, and
    % reset the end point azimuth to match azLimit exactly.
    azLimit = az(1) + [0 2*pi];
    az(end) = azLimit(2);
else
    % Otherwise, try to avoid breaking the ring into parts.
    azLimit = min(az) + [0 2*pi];
end

% Treat the ring as a polyline, and trim in azimuth (which is cyclical in
% the same way that longitude is). Rather than remove segments of the ring,
% this operation will break it into parts as needed, with each part wrapped
% into the interval defined by azLimit and terminated at the limits of that
% interval.
[rng, az] = trimPolylineToLonlim(rng, az, azLimit);

% Clean out any vertical strips left over from polyline trimming.
[rng, az] = removeResidualsAtLimits(rng, az, azLimit);

% As needed, reconnect dangling endpoints along the circumference of a
% rectangle with xLimit equal to azLimit and yLimit (rngLimit) values
% chosen to comfortable enclose the entire ring. The first and last values
% of az already match azLimit exactly, so set tolSnap to zero to suppress
% further snapping.
rngLimit = [0 2*max(rng)];
tolSnap = 0;
[az, rng] = closePolygonInRectangle(az, rng, azLimit, rngLimit, tolSnap);

% Wrap the test point azimuths into the half-open interval
% [azLimit(1) azLimit(2)).
azTest = azLimit(1) + mod(azTest - azLimit(1), 2*pi);

% Ready now to use the MATLAB inpolygon function.
tfLeft = inpolygon(azTest, rngTest, az, rng);

% To be certain, repeat using the half-open interval (azLimit(1) azLimit(2)].
azTest(azTest == azLimit(1)) = azLimit(2);
tfRight = inpolygon(azTest, rngTest, az, rng);

% Combine results.
tf = tfLeft | tfRight;

%-----------------------------------------------------------------------

function [rng, az] = removeResidualsAtLimits(rng, az, azLimit)
% If the process of trimming a ring in azimuth has left any strips of
% constant range along either limit, remove such strips. Adapted from the
% subfunction of the same name in private/trimPolygonToQuadrangle.

[first, last] = internal.map.findFirstLastNonNan(rng);
iRemove = false(size(rng));
for k = 1:numel(first)
    if all(az(first(k):last(k)) <= azLimit(1)) || ...
       all(az(first(k):last(k)) >= azLimit(2))
        iRemove(first(k):last(k)) = true;
    end
end
rng(iRemove) = [];
az(iRemove) = [];

%--------------------------------------------------------------------------

function [rng, az] = latlon2rngaz(lat, lon, latCenter, lonCenter)

% Transforms latitude-longitude to range-azimuth and unwrap the
% azimuth angles.  Note special handling when the center is a pole.
if latCenter >= pi/2
    rng = pi/2 - lat;
    az = -lon;
elseif latCenter <= -pi/2
    rng = lat + pi/2;
    az  = lon;
else
    [rng, az] = distance(latCenter, lonCenter, lat, lon, 'radians');
end
az = unwrapMultipart(az);

% For rings in which azimuths don't wrap, restore exact match between
% first and last elements.
[first,last] = internal.map.findFirstLastNonNan(az);
nonWrapping = abs(az(last) - az(first)) < pi;
az(last(nonWrapping)) = az(first(nonWrapping));
rng(last(nonWrapping)) = rng(first(nonWrapping));

[az,rng] = shiftInnerRingsBy2PiN(az,rng);

%--------------------------------------------------------------------------

function [lat, lon] = rngaz2latlon(rng, az, latCenter, lonCenter)

% Transforms range-azimuth to latitude-longitude.
if latCenter >= pi/2
    lat = pi/2 - rng;
    lon = -az;
elseif latCenter <= -pi/2
    lat = rng - pi/2;
    lon = az;
else
    [lat, lon] = reckon(latCenter, lonCenter, rng, az, 'radians');
end
lon = wrapToPi(lon);

%--------------------------------------------------------------------------

function [az,rng] = shiftInnerRingsBy2PiN(az,rng)
% (AZ,RNG) is an azimuth-range representation of a set of polygon parts,
% or rings. AZ is in radians. RNG is non-negative. Both inputs are
% NaN-separated and NaN-terminated, and have identical sizes and NaN
% locations. Neither is empty. They have already been unwrapped in
% azimuth.
%
% The rings represented by (AZ,RNG) may be unwrapped in azimuth
% independently of each other. This means that a counter-clockwise ring
% that is actually an inner ring of a different, clockwise ring, may be
% shifted in azimuth relative to its outer ring, by an amount equal to
% 2*Pi*N, where N is nonzero integer.  (We expect N to be 1 or -1, but
% also check for 2 and -2).
%
% This function identifies such situations and undoes the shifts. In its
% output every counter-clockwise ring is contained within a clockwise
% ring, when viewed in a right-handed pseudo-Cartesian system in which
% AZ is equated with X and RNG is equated with Y, as long as there is a
% shift that makes this possible. Given otherwise valid inputs, this
% function should not fail unless the topology is corrupt.

% Find the first and last vertices for each ring.  These will help us
% reference and manipulate the coordinates using only "in-place"
% operations.
[first, last] = internal.map.findFirstLastNonNan(az);

% Try 5 possible shifts.
shifts = 2*pi * [-2 -1 0 1 2];

% Determine which rings are cw and which are ccw
cw = ispolycw(az,rng);
cwIndex  = find(cw);
ccwIndex = find(~cw);
M = numel(ccwIndex);

% Construct M-by-5 arrays of points to check with INPOLYGON:
% One row per ccw ring and one column per possible shift.
checkAz  = az( first(~cw),[1 1 1 1 1]) + shifts(ones(M,1),:);
checkRng = rng(first(~cw),[1 1 1 1 1]);

% Keep track of the cw polygon with which each ccw polygon has been
% associated.
containingPoly = NaN + zeros(M,1);

% Iterate over the set of cw polygons
for p = 1:numel(cwIndex)
    % Check all the ccw rings for inclusion in the current
    % cw ring, given shifts of N = -2, -1, 0, 1, and 2.
    k = cwIndex(p);
    in = inpolygon(checkAz, checkRng, ...
                   az(first(k):last(k)), rng(first(k):last(k)));

    % Loop over the ccw rings contained in the current cw ring.
    indx = find(any(in,2));
    for q = 1:numel(indx)
        j = indx(q);               
        if cwPolygonIsNested(k,containingPoly(j),az,rng,first,last)               
            shift = shifts(in(j,:));
        
            % There should be exactly one 2*pi*N shift that puts
            % the current ccw ring inside the current cw ring.
            assert(isscalar(shift), ['map:' mfilename ':internalError'], ...
                'Mapping Toolbox internal error: Expected %s to be a scalar.', ...
                'SHIFT')
            
            m = ccwIndex(j);
            az(first(m):last(m)) = az(first(m):last(m)) + shift;
            containingPoly(j) = k;
        end
    end
end

%-----------------------------------------------------------------------

function tf = cwPolygonIsNested(k,m,az,rng,first,last)
% True if the cw polygon with index k is nested within the cw polygon
% with index m, or if we haven't yet identified a containing polygon --
% in which case m would equal NaN.

if isnan(m)
    tf = true;
else
    tf = inpolygon(az(first(k)),         rng(first(k)), ...
                   az(first(m):last(m)), rng(first(m):last(m)));
end

%-----------------------------------------------------------------------

function A = areaPolar(th, r)
% Area between curve and origin in polar coordinate system
%
%   Given a set of ordered vertices in polar coordinates defined by
%   vectors th and r, compute the area between the origin and the curve
%   defined by interpolating linearly for r in terms of th between each
%   adjacent pair of vertices.
%
%   Example
%   -------
%   r =  [1   1     1      2     2      2     2     2    1 ];
%   th = [0 pi/4 3*pi/4 3*pi/4   pi  5*pi/4 3*pi/2  0    0 ];
%   polar(th,r,'ro')
%   A = areaPolar(th,r)
%   A_expected = (3/8)*pi + (5/8)*(4*pi)

% For each adjacent pair of points, (r1,th1) and (r2, th2), define the
% curve that connects them via linear interpolation of r in terms of th
% as follows:
%
%     r(th) = (r2 * (th - th1) + r1 * (th2 - th)) / (th2 - th1).
%
% Integrate to find the area of this curve:
%
%     A12 = integral of ((r(th))^2)/2 d(th) from th = th1 to th = th2.
%
% The value of this definite integral can be expressed analytically as:
%
%     A12 = (th2 - th1) * (r1^2 + r1*r2 + r2^2) / 6.
%
% The preceding integration was performed using the following script
% (which requires Symbolic Toolbox to run):
% 
%     syms r1 r2 th1 th2 r th
%     r = (r2*(th - th1) + r1*(th2 - th)) / (th2 - th1);
%     int(r^2/2, th, th1, th2)

r1 = r(1:end-1);
r2 = r(2:end);
mean_r_squared = (r1.^2 + r1.*r2 + r2.^2)/3;
delta_th = diff(unwrap(th));
A = sum(mean_r_squared .* delta_th) / 2;
