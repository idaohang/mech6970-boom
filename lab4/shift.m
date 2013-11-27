function y = shift(x,n)

% Move all vector values up towards the front by n indices
% 
% USAGE:
%   [2 3 4 1] = shift([1 2 3 4], 1)
%   [4 1 2 3] = shift([1 2 3 4], -1)
% 
% @author Robert Cofield
% 

if isnan(n)
  warning('Received NaN shift value. Returning input as output.');
  y = x;
  return
elseif round(n) ~= n
  error('Non-integer shift value received');
end

if abs(n) > length(x)
  n = sign(n)*mod(abs(n),length(x));
end

if n == 0
%   warning('Received `n` value of 0. Returning input vector as output');
  y = x;
  return
elseif n > 0
  y = [x(1+n:end) x(1:n)];
elseif n < 0
  n = -n;
  y =  [x(end-n+1:end) x(1:end-n)];
end

end