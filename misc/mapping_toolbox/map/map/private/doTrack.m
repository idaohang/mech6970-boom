function [outPhi, outLambda] = doTrack(trackstr, phi, lambda ,ellipsoid, npts)
% Core computations performed by TRACK.  All angles are in units of
% radians.  TRACKSTR can be 'gc' or 'rh'.

% Copyright 2006 The MathWorks, Inc.
% $Revision: 1.1.6.1 $  $Date: 2006/12/10 20:04:42 $

%  Ensure that phi and lambda are column vectors
phi = phi(:);
lambda = lambda(:);

%  Compute vectors of start and end points
startlats = phi(1:length(phi)-1);
endlats   = phi(2:length(phi));
startlons = lambda(1:length(lambda)-1);
endlons   = lambda(2:length(lambda));

[outPhi,outLambda] = doTrack2(...
    trackstr, startlats, startlons, endlats, endlons, ellipsoid, npts);

%  Link all tracks into a single NaN clipped vector
[r,c] = size(outPhi); %#ok<NASGU>
outPhi(r+1,:) = NaN;
outLambda(r+1,:) = NaN;
outPhi = outPhi(:);
outLambda = outLambda(:);
