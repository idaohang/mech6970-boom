function [Z, refvec] = spzerom(latlim, lonlim, scale)
%SPZEROM  Construct sparse regular data grid of 0s
%
%   [Z, REFVEC] = SPZEROM(LATLIM, LONLIM, SCALE) constructs a sparse
%   regular data grid consisting entirely of 0s.  The two-element
%   vectors LATLIM and LONLIM define the latitude and longitude limits
%   of the grid, in degrees.  They should be of the form [south north]
%   and [west east], respectively.  The number of rows and columns per
%   degree is set by the scalar value SCALE.  REFVEC is the
%   three-element referencing vector for the data grid.
%
%   See also NANM, ONEM, ZEROM, SPARSE.

% Copyright 1996-2006 The MathWorks, Inc.
% $Revision: 1.10.4.3 $  $Date: 2006/05/24 03:36:44 $

[nrows, ncols, refvec] = sizem(latlim, lonlim, scale);
Z = sparse(nrows, ncols);
