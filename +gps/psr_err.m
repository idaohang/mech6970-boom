function varargout = psr_err(psr, sv_pos, user_pos, varargin)
% PSR_ERR calculates the error in pseudorange measurements given and
% estimate pseudoranges and estimates for user position, SV position, and
% user clock bias
% 
% USAGE:
%   e = psr_err(psr, sv_pos, user_pos, clock_bias)
%   [e,psr_est] = psr_err(psr, sv_pos, user_pos, clock_bias)
% 
% INPUTS:
%   psr: 1xn vector of pseudoranges (meters)
%   sv_pos: nx3 matrix of SV positions (ECEF meters) --- x,y,z
%   user_pos: 1x3 vector for ECEF user position (meters)
%   clock_bias: estimate for user's
% 
% OUTPUTS:
%   e: 1xn vector of psr corrections
%     e = psr - psr_est
%   psr_est: 1xn vector of psr estimates

clock_bias = varargin{1};

% number of SV's
nsv = length(psr);

psr_est = sqrt(sum((sv_pos-repmat(user_pos,nsv,1)).^2, 2)) + c*clock_bias;
e = psr - psr_est;

varargout{1} = e;
if nargout == 2
  varargout{2} = psr_est;
end

end