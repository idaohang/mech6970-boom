function checkboundingbox(...
               bbox, function_name, variable_name, argument_position)
%CHECKBOUNDINGBOX Check validity of bounding box array.
%   CHECKBOUNDINGBOX(...
%              BBOX, FUNCTION_NAME, VARIABLE_NAME, ARGUMENT_POSITION)
%   ensures that the bounding box array is a 2-by-2 array of double with
%   real, finite values, and that in each column the second value always
%   exceeds the first.

% Copyright 1996-2007 The MathWorks, Inc.
% $Revision: 1.1.6.2 $  $Date: 2007/12/10 21:38:56 $

component = 'map';

checkinput(bbox, {'double'},{'real','nonnan'},...
           function_name,variable_name,argument_position);
        
if size(bbox,1) ~= 2 || size(bbox,2) ~= 2 || ndims(bbox) ~= 2
    eid = sprintf('%s:%s:invalidBBoxSize',component,function_name);
    msg1 = sprintf('Function %s expected its %s input argument, %s,', ...
                   upper(function_name), num2ordinal(argument_position), ...
                   variable_name);
    msg2 = sprintf('to be a 2-by-2 array.');
    error(eid, '%s\n%s', msg1, msg2);
end

if ~all(bbox(1,:) <= bbox(2,:))
    eid = sprintf('%s:%s:invalidBBoxOrder',component,function_name);
    msg1 = sprintf('Function %s expected its %s input argument, %s,', ...
                   upper(function_name), num2ordinal(argument_position), ...
                   variable_name);
    msg2 = sprintf('to have element (2,k) greater than or equal to element (1,k).');
    error(eid, '%s\n%s', msg1, msg2);
end
