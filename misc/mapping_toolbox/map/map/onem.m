function [Z, refvec] = onem(latlim, lonlim, scale)
%ONEM  Construct regular data grid of 1s
%
%   [Z, REFVEC] = ONEM(LATLIM, LONLIM, SCALE) constructs a regular
%   data grid consisting entirely of 1s.  The two-element vectors LATLIM
%   and LONLIM define the latitude and longitude limits of the grid, in
%   degrees.  They should be of the form [south north] and [west east],
%   respectively.  The number of rows and columns per degree is set by
%   the scalar value SCALE.  REFVEC is the three-element referencing
%   vector for the data grid.
%
%   See also NANM, SPZEROM, ONES, ZEROM.

% Copyright 1996-2006 The MathWorks, Inc.
% $Revision: 1.10.4.3 $  $Date: 2006/05/24 03:35:14 $

[nrows, ncols, refvec] = sizem(latlim, lonlim, scale);
Z = ones(nrows, ncols);
