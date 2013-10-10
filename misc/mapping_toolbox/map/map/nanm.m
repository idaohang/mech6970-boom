function [Z, refvec] = nanm(latlim, lonlim, scale)
%NANM  Construct regular data grid of NaNs
%
%   [Z, REFVEC] = NANM(LATLIM, LONLIM, SCALE) constructs a regular
%   data grid consisting entirely of NaNs.  The two-element vectors
%   LATLIM and LONLIM define the latitude and longitude limits of the
%   grid, in degrees.  They should be of the form [south north] and
%   [west east], respectively.  The number of rows and columns per
%   degree is set by the scalar value SCALE.  REFVEC is the
%   three-element referencing vector for the data grid.
%
%   See also SPZEROM, ONEM, ZEROM.

% Copyright 1996-2006 The MathWorks, Inc.
% $Revision: 1.10.4.3 $  $Date: 2006/05/24 03:35:05 $

[Z, refvec] = zerom(latlim, lonlim, scale);
Z = NaN + Z;
