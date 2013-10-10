function RGB = ind2rgb8(X, CMAP)
%IND2RGB8 Convert indexed image to uint8 RGB image
%
%   RGB = IND2RGB8(X,CMAP) creates an RGB image of class uint8.  X must be
%   uint8, uint16, or double, and CMAP must be a valid MATLAB colormap.
%
%   Example 
%   -------
%      % Convert the 'concord_ortho_e.tif' image to RGB.
%      [X, cmap] = imread('concord_ortho_e.tif');
%      RGB = ind2rgb8(X, cmap);
%      R = worldfileread('concord_ortho_e.tfw');
%      mapshow(RGB, R);
%
%   See also IND2RGB.

% Copyright 1996-2007 The MathWorks, Inc.
% $Revision: 1.1.6.3 $  $Date: 2007/06/04 21:11:51 $ 

RGB = ind2rgb8c(X, CMAP);
