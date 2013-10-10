function str = leadblnk(str,char)
%LEADBLNK  Deletes leading characters common to all rows of a string matrix
%
%  s = LEADBLNK(s0) deletes the leading spaces common to all
%  rows of a string matrix.
%
%  s = LEADBLNK(s0,char) deletes the leading characters corresponding
%  to the scalar char.
%
%  See also  NUM2STR, SHIFTSPC.

% Copyright 1996-2007 The MathWorks, Inc.
% $Revision: 1.9.4.2 $  $Date: 2007/11/09 20:27:37 $
% Written by:  E. Byrns, E. Brown

error(nargchk(1, 2, nargin, 'struct'))

if nargin == 1
    char = [];
end

%  Empty argument tests
if isempty(str)
    return
end
if isempty(char)
    char = ' ';
end

%  Eliminate null characters
indx = find(str==0);
if ~isempty(indx)
    str(indx) = char;
end

%  Eliminate leading spaces
if size(str,1) ~= 1;
    spaces = sum(str ~= char);    %  Not spaces in each row
else
    spaces = (str ~= char);       %  Single string entry
end

indx = find(spaces~=0, 1);           %  0 sum column-wise means char is in every row
if isempty(indx)
    indx = length(spaces)+1;
end    %  All char in matrix
str(:,1:(indx-1)) = [];
