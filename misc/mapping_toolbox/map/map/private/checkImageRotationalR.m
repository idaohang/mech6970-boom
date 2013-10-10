function [dataArgs, numArgs] = ...
    checkImageRotationalR(isImage, dataArgs, numArgs)
%CHECKIMAGEROTATIONALR Update dataArgs if R is not rectilinear
%
%   [dataArgs, numArgs] = checkImageRotationalR(isImage, dataArgs, numArgs)
%   checks if dataArgs contains an image and a referencing matrix that is 
%   not rectilinear. dataArgs must have one of the following forms:
%      {}, {I}, {I, R}, {I, CMAP, R}, {C1, C2, I}, or {C1, C2, I, CMAP}.
%   dataArgs is changed only under the following conditions: I is an image
%   and R is a non-rectilinear reference matrix. In this case, {I, R} is
%   replaced with {C1, C2, I} and {I, CMAP, R} is replaced with
%   {C1, C2, I, CMAP}. A referencing matrix R is non-rectilinear if it has
%   any non-zero diagonal elements.

% Copyright 2010 The MathWorks, Inc.
% $Revision: 1.1.6.1 $  $Date: 2010/03/22 03:51:56 $

twoOrThreeArgs = numArgs == 2 || numArgs == 3;
isRefMat = isequal(size(dataArgs{numArgs}),[3 2]);
if isImage && twoOrThreeArgs && isRefMat && any(diag(dataArgs{numArgs}))
    
    % Syntax: (I, R) or (I, CMAP, R)
    % R is rotational, create new dataArgs
    [C1, C2] = pixcenters(dataArgs{numArgs}, size(dataArgs{1}));
    
    % Remove R
    dataArgs(numArgs) = [];
    
    % New dataArgs: (C1, C2, I, ...)
    dataArgs = [{C1}, {C2}, dataArgs];
    numArgs = numArgs + 1;
end
