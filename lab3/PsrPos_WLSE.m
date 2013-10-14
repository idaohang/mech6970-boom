function [pos_soln, bias_soln, itx] = PsrPos_LSE(psr, svpos, R, P0, pos0, bias0, pos_tol, maxit)

% Estimate (Recursive Least Squares) ECEF position using pseudorange estimates
% 
% USAGE:
%   [pos_soln, bias_soln, itx] = PsrPos_LSE(psr, svpos, R, P0, pos0, bias0, pos_tol, maxit)
%   
% INPUTS:
%   psr: Nx1 matrix of Pseudoranges from each satellite
%   svpos: Nx3 matrix of ECEF SV Positions
%   R: measurement covariance matrix
%   P0: Initial Estimate error covariance matrix
%   pos0: Initial position guess
%   bias0: Intial clock bias guess
%   pos_tol: iterate when norm of position change is below this
%   maxit: maximum number of iterations
% 
% OUTPUTS:
%   pos_soln: 3x1 ECEF position solution
%   bias_soln: user clock bias solution (m)
% 

P = P0;
est = [pos0;bias0];

mv = pos_tol+1;
ct = 1;

while mv > pos_tol
  est_old = est;
  x0 = est_old(1:3); % set new initial guess to last solution
  b0 = est_old(4);
  psr0 = zeros(length(svpos),1);
  % calculate psr guess
  for k = 1:length(svpos)
    psr0(k) = norm(svpos(k,:)-x0', 2) + b0;
  end
  % measurement inputs
  dpsr = psr' - psr0;  
  % calculate design matrix
  G = [gps.calc_geometry_matrix(est(1:3), svpos) ones(length(svpos),1)];  
  % calculate gain
  K = P*G'*inv(G*P*G'+R);  
  % new estimate of offset from guess
  %   dx = dx + K*(dpsr - G*dx);  
  dx = K*dpsr;
  % resultant position and bias estimates
  est = [x0;b0] + dx;  
  % New estimation error covariance
  P = (eye(4) - K*G)*P;    
  mv = norm(dx(1:3), 2);
  if ct>maxit
    break;
  end  
  ct = ct+1;
end

% outputs

pos_soln = est(1:3);
bias_soln = est(4);
if nargout > 2
  itx = ct-1;
end

