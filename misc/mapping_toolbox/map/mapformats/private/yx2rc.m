function [row,col] = yx2rc(y,x,y1,x1,yperrow,xpercol)

%YX2RC X and Y coordinates to row and column indices
% [rowgrat,colgrat] = YX2RC(ygrat,xgrat,y1,x1,yperrow,xpercol) computes row 
% and column indices from y and x coordinates. ygrat and xgrat are the 
% cartesian coordinates for a grid with the upper left corner at coordinates 
% x1 and y1. yperrow and xpercol are the y and x increment per row and column. 
% The outputs rowgrat and colgrat correspond to matrix locations ygrat and 
% xgrat.
%
% See also RC2YX, READMTX.

% Copyright 1996-2007 The MathWorks, Inc.
% $Revision: 1.1.6.2 $  $Date: 2007/10/10 20:49:44 $

if nargin~=6
   error(['map:' mfilename ':mapformatsError'], 'Incorrect number of arguments')
end

row = ceil( 1 + (y - y1)/yperrow );
col = ceil( 1 + (x - x1)/xpercol );




