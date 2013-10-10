function RGB = checkImage(mapfilename, A, cmap, imagePos, cmapPos)
%CHECKIMAGE Check image inputs
%
%   RGB = checkImage(MAPFILENAME, A, CMAP, imagePos, cmapPos) validates the
%   image input, A. If A is not an RGB image, then A is converted to an RGB
%   image by using IND2RGB8 if CMAP is not empty or by replicating the 
%   matrix. If A is logical it is converted to uint8. imagePos and cmapPos
%   are the position numbers in the command line for the arguments.

% Copyright 2010 The MathWorks, Inc.
% $Revision: 1.1.6.3 $  $Date: 2010/09/24 14:33:45 $

if islogical(A)
   u = uint8(A);
   u(A) = 255;
   A = u;
end

if ~isempty(cmap)
    internal.map.checkcmap(cmap, mapfilename, 'CMAP', cmapPos);
    attributes = {'2d', 'real', 'nonsparse', 'nonempty'};
    checkinput(A, {'numeric'}, attributes, mapfilename, 'I or X or RGB', ...
        imagePos);
    A = ind2rgb8(A, cmap);
else
    attributes = {'real', 'nonsparse', 'nonempty'};
    checkinput(A, {'numeric'}, attributes, mapfilename, 'I or X or RGB', ...
        imagePos);
    
    if ndims(A) == 2
        A = repmat(A,[1 1 3]);
    elseif ndims(A) ~= 3
        eid = sprintf('%s:%s:invalidImageDimension', getcomp, mapfilename);
        msg = sprintf('%s','Image dimension must be 2 or 3.');
        error(eid, '%s',msg)
    end   
end
RGB = checkRGBImage(A);

%--------------------------------------------------------------------------

function RGB = checkRGBImage(RGB)

% RGB images can be only uint8, uint16, or double
if ~isa(RGB, 'double') && ...
      ~isa(RGB, 'uint8')  && ...
      ~isa(RGB, 'uint16')
   eid = sprintf('%s:%s:invalidRGBClass', getcomp, mfilename);
   error(eid, 'RGB images must be uint8, uint16, or double.');
end

if size(RGB,3) ~= 3
   eid = sprintf('%s:%s:invalidRGBSize', getcomp, mfilename);
   error(eid, 'RGB images must be size M-by-N-by-3.');
end

% Clip double RGB images to [0 1] range
if isa(RGB, 'double')
   RGB(RGB < 0) = 0;
   RGB(RGB > 1) = 1;
end
