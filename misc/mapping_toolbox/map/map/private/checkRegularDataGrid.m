function [Z, R] = checkRegularDataGrid(Z, R, fcnName)
%CHECKREGULARDATAGRID Check regular data grid inputs
%
%   [Z, R] = checkRegularDataGrid(Z, R, fcnName) validates the regular data
%   grid, defined by the inputs, Z and R. A referencing matrix is returned
%   in R. fcnName is used for constructing an error message. 

% Copyright 2010 The MathWorks, Inc.
% $Revision: 1.1.6.1 $  $Date: 2010/03/22 03:51:58 $

classes = {'double', 'single'};
attributes = {'real', '2d', 'nonempty'};
checkinput(Z, classes, attributes, fcnName, 'Z' ,1);
R = checkRefObj(fcnName, R, size(Z), 2);