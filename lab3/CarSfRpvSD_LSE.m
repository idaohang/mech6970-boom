function [rpv_soln, bias_soln, itx, dop] = CarSfRpvSD_LSE(carsd, afixed, wavelength, svpos, rpv0, bias0, pos_tol, maxit)

% Estimate (Recursive Least Squares) ECEF position using single differenced
% carrier phase along with (presumably) fixed integer ambiguity estimates.
% Single frequency carrier data is supported.
% 
% USAGE:
%   [rpv_soln, bias_soln, itx, dop] = CarSfRpvSD_LSE(carsd, afixed, wavelength, svpos, pos0, bias0, pos_tol, maxit)
%   
% INPUTS:
%   carsd: Nx1 matrix of single-differenced carrier phase ranges
%   afixed: Nx1 vector of carrier phase ambiguity estimates
%   wavelength: length of the carrier waveform
%   svpos: Nx3 matrix of ECEF SV Positions
%   rpv0: Initial relative position guess
%   bias0: Intial clock bias guess
%   pos_tol: iterate until norm of relative position change is below this
%   maxit: maximum number of iterations
% 
% OUTPUTS:
%   rpv_soln: 3x1 ECEF relative position vector solution
%   bias_soln: user clock bias difference solution (m)
%   dop: matrix of various dilutions of precision
%     [HDOP, VDOP, PDOP, TDOP, GDOP]
% 

est = [rpv0;bias0];

mv = pos_tol+1;
ct = 1;

% The original measurement
y_orig = carsd - wavelength*afixed;

% calculate design matrix
G = [calc_geometry_matrix(est(1:3), svpos) ones(length(svpos),1)];  

while mv > pos_tol
  
  est_old = est;
  x0 = est_old(1:3); % set new initial guess to last solution
  b0 = est_old(4);
  
  % calculate the measurement that would result from this being the state
  y0 = G*est_old;  
  
  % measurement inputs
  dy = y_orig - y0;  
  
  % calculate design matrix
  G = [calc_geometry_matrix(est(1:3), svpos) ones(length(svpos),1)];  
  % new estimate of offset from guess
  dx = pinv(G)*dy;
  % resultant position and bias estimates
  est = [x0;b0] + dx;  

  mv = norm(dx, 2);
  if ct>maxit
    break;
  end  
  ct = ct+1;
end

% outputs

rpv_soln = est(1:3);
bias_soln = est(4);
if nargout > 2
  itx = ct-1;
  if nargout > 3
    dop = calcDOP(G);
  end
end

