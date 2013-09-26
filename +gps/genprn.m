function prn = genprn(n, nchips, scale, freq)
% prn = genprn(n, len, scale, freq)
% 
% GENPRN generate a pseudorandom noise (PRN) sequence, otherwise known as
% C/A code, or Gold code.
% 
% INPUTS:
%   n: PRN #'s (1-37 supported)
%   len: desired length of the PRN 
%   scale: 2 element vector of desired values for low and high bits,
%     respectively
%   freq: upsamples by multiplying the chiprate
% 
% OUTPUTS:
%   prn: length(n) x len matrix. Each SV's PRN is contained in the rows.
% 
% USAGE:
%   prn = gps.genprn([12 19], 1023)
%   prn = gps.genprn([12 19], 1023, [1 -1])
%   prn = gps.genprn([12 19], 1023, [-1 1], 3)
% 

% inputs
if nargin < 4
  freq = false;
  if nargin < 3
    scale = false;
  end
end

% TODO input error checking here

% load PRN tap indices
gps.constants
% number of satellites for which the prn will be calculated
nsv = length(n);
% Registers
G1 = ones(1,10);
G2 = ones(1,10);
% which indices to use for the G2 register
S = tap(n,:);
% output
prn = zeros(nsv,nchips);

for k = 1:nchips
  % compute output of G2 register
  G2k = xor(G2(S(:,1)),G2(S(:,2)));
  % compute next PRN bit
  prn(:,k) = xor(G1(end),G2k);
  % Compute new values for registers
  G1 = [ xor(G1(3),G1(10)) , G1(1:end-1) ];
  G2 = [ mod(sum(G2([2,3,6,8,9,10]),2),2) , G2(1:end-1) ];
end

% Scaling
if scale
  prn_ = prn;
  prn(prn_==0) = scale(1);
  prn(prn_==1) = scale(2);
end

% Upsampling
if freq
  prn_ = prn;
  prn = zeros(1,freq*nchips);
  for k = 1:nchips
    prn((k-1)*freq+1:(k-1)*freq+freq) = prn_(k);
  end
end

end