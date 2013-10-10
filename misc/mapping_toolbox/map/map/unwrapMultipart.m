function p = unwrapMultipart(p)
%unwrapMultipart Unwrap vector of angles with NaN-delimited parts
%
%   UNWRAPPED = unwrapMultipart(P) unwraps a row or column vector of
%   azimuths, longitudes, or phase angles.  If P is separated into
%   multiple parts delimited by values of NaN, then each part is
%   unwrapped independently.  If P has only one part, then the result is 
%   equivalent to UNWRAP(P).  Input and output are both in radians. The
%   output is the same size as the input and has NaNs in the same
%   locations.
%
%   See also UNWRAP.

% Copyright 2007-2009 The MathWorks, Inc.
% $Revision: 1.1.6.1 $  $Date: 2009/03/09 19:16:24 $

% The following implementation was generalized from subfunction
% LocalUnwrap in wrap.m

% Set cutoff/tolerance to the default used by unwrap.
cutoff = pi;

% Remember input size and shape.
sizep = size(p);

% Make sure that p is a column vector.
p = p(:);

% Incremental phase variations.
dp = [0; diff(p,1,1)];

% Incremental phase corrections.
dp_corr = wrapToPi(dp) - dp;

% Zero correction when incremental variation is less than tolerance.
dp_corr(abs(dp) < cutoff) = 0;

% Zero out NaNs in dp_corr so that they don't block the cumulative sum.
dp_corr(isnan(dp_corr)) = 0; 

% Find NaN positions in p, and keep only those that are either
% (a) isolated or (b) the final element in a contiguous sequence of NaN.
n = find(isnan(p));
n(diff(n) == 1) = [];

% If NaNs are present, adjust corrections so that their cumulative sum
% starts over at zero following each sequence of contiguous NaNs.
if ~isempty(n)
    rawcumsum = cumsum(dp_corr,1);
    dp_corr(n) = dp_corr(n) - diff([0; rawcumsum(n)]);
end

% Integrate and apply the adjusted corrections.
p = p + cumsum(dp_corr,1);

% Restore input shape.
p = reshape(p,sizep);
