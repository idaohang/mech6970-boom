function [ygrat,xgrat] = rc2yx(rowgrat,colgrat,y1,x1,yperrow,xpercol)

%RC2YX Row and column indices to x and y coordinates
%
% [ygrat,xgrat] = RC2YX(rowgrat,colgrat,y1,x1,yperrow,xpercol) computes y and x
%  coordinates from row and column indices. rowgrat and colgrat are matrices 
%  of the same size containing integer row and column indices for a grid with 
% the upper left corner at the Cartesian coordinates x1 and y1. yperrow and 
% xpercol are the y and x increment per row and column. The outputs ygrat and 
% xgrat are the cartesian coordinates corresponding to matrix locations rowgrat
% and colgrat.
%
% See also YX2RC, READMTX.

% Copyright 1996-2007 The MathWorks, Inc.
% $Revision: 1.1.6.2 $  $Date: 2007/10/10 20:49:39 $

if nargin~=6
	error(['map:' mfilename ':mapformatsError'], 'Incorrect number of arguments')
end


ygrat = (rowgrat-1)*yperrow + y1;
xgrat = (colgrat-1)*xpercol + x1;

