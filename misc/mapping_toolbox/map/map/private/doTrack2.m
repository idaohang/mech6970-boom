function [phiTrack, lambdaTrack] ...
    = doTrack2(str, phi1, lambda1, phi2, lambda2, ellipsoid, npts)
% Core computations performed by TRACK2.  All angles are in radians.

% Copyright 2006 The MathWorks, Inc.
% $Revision: 1.1.6.1 $  $Date: 2006/12/10 20:04:43 $

% Ensure that inputs are column vectors
phi1 = phi1(:);
phi2 = phi2(:);
lambda1 = lambda1(:);
lambda2 = lambda2(:);

% Set tolerance
epsilon = epsm('radians');

%  If a track starts at a pole, ensure it traverses
%  straight down the appropriate longitude
startingAtPole = (abs(phi1) >= pi/2-epsilon);
lambda1(startingAtPole) = lambda2(startingAtPole);

%  DISTANCE and RECKON handle ellipsoid input of [0 0],
%  so we need to take care of that here, also.
if ellipsoid(1) == 0
    ellipsoid(1) = 1;
end

%  Compute distance and course
isSphere = (numel(ellipsoid) < 2 || ellipsoid(2) == 0);
if isSphere && strcmp(str,'gc')
    [rng, az] = greatcircleinv(phi1, lambda1, phi2, lambda2, ellipsoid(1));
else
    [rng, az] = distance(str, phi1, lambda1, phi2, lambda2, ellipsoid, 'radians');
end

%  Compute the tracks -- to avoid unit conversions and other
%  inefficiencies, perform just the essential track1 computation,
%  producing results equivalent to the following call:
%
% [phiTrack,lambdaTrack] ...
%     = track1(str, phi1, lambda1, az, rng, ellipsoid, 'radians', npts);

%  Use real(npts) to avoid a cumbersome warning for complex n in linspace
npts = real(npts);
bigsize = [size(rng,1), npts];
bigrng = zeros(bigsize);
for i = 1:size(rng,1)
	bigrng(i,:) = linspace(0, rng(i), npts);
end

%  Compute the tracks
%  Each track occupies a row of the output matrices.
bigPhi    = phi1(    :, ones([1,npts]));
bigLambda = lambda1( :, ones([1,npts]));
bigaz     = az(      :, ones([1,npts]));

% Expand the following call to reckon to by-pass argument checking and
% angle units conversion:
%
% [phiTrack, lambdaTrack] ...
%     = reckon(str, bigPhi(:), bigLambda(:), bigrng(:), bigaz(:), ellipsoid, 'radians');

if strcmp(str,'gc')
    if ellipsoid(2) ~= 0
        [phiTrack, lambdaTrack] = geodesicfwd(...
            bigPhi(:), bigLambda(:), bigaz(:), bigrng(:), ellipsoid);
    else
        [phiTrack, lambdaTrack] = greatcirclefwd(...
            bigPhi(:), bigLambda(:), bigaz(:), bigrng(:), ellipsoid(1));
    end
elseif strcmp(str,'rh')
    if ellipsoid(2) ~= 0
        [phiTrack, lambdaTrack] = rhumblinefwd(...
            bigPhi(:), bigLambda(:), bigaz(:), bigrng(:), ellipsoid);
    else
        [phiTrack, lambdaTrack] = rhumblinefwd(...
            bigPhi(:), bigLambda(:), bigaz(:), bigrng(:), ellipsoid(1));
    end
else
    eid = sprintf('%s:%s:invalidTrackString',getcomp,mfilename);
    error(eid,'Invalid track string: %s',str)
end

lambdaTrack = npi2pi(lambdaTrack,'radians');

%  Restore shape, then transpose the reckon results so that each track
%  occupies one column of the output matrices.
phiTrack = reshape(phiTrack,bigsize)';
lambdaTrack = reshape(lambdaTrack,bigsize)';
