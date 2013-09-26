function [seq, record] = lsr(seed, n)
% LSR performs the action of a variable-bit long shift register, as done to 
% produce C/A code. It uses modulo 2 addition.
% Misra,Enge pg 62.
% 
% INPUTS:
%   seed: the initial sequence
%   n: number of epochs (iterations) to perform the shift operation. T
% 
% OUTPUTS:
%   record: initial seed with subsequent rows containing intermediate
%     values and what would have been the next bitstring
%     [
%   seq: resulting final sequence
% 
% USAGE:
%   [seq, record] = lsr(seed, n)
% 
% EXAMPLE:
%   seed = [1 1 0 0 1]
%   len
% 
if (n<1) || (n~=floor(n))
  error('Improper epoch number `n`');
end

len = length(seed); % how many bits for each entry in record
record = zeros(n,len);
record(1,:) = seed;
seq = ones(1,n);

for k = 1:n
  seq(k) = record(k,end);
  record(k+1,2:end) = record(k,1:end-1);
  record(k+1,1) = xor(record(k,2),record(k,end));
end



end