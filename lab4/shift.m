function y = shift(x,n)

% Move all vector values up towards the front by n indices
% 
% USAGE:
%   [2 3 4 1] = shift([1 2 3 4], 1)
%   [4 1 2 3] = shift([1 2 3 4], -1)
% 
% @author Robert Cofield
% 

if round(n) ~= n
  error('Non-integer shift value received');
end

if abs(n) > length(x)
  n = sign(n)*mod(abs(n),length(x));
end

if n == 0 
  y = x;
elseif n > 0
  y = [x(1+n:end) x(1:n)];
elseif n < 0
  n = -n;
  y =  [x(end-n+1:end) x(1:end-n)];
end

end