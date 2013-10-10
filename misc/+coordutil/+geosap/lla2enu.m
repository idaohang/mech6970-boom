function [enu] = lla2enu(lla,origin_lla)
%LLA2ENU Transforms Lat/Lon/Alt to East/North/Up.
%
%   [enu] = lla2enu(lla,origin_lla)
%
%   Inputs:
%   lla = n x 3 matrix of lat/lon/alt data to be transformed
%   origin_lla = 1 x 3 array describing the desired origin of the ENU coordinate frame
%       both inputs have units [radians radians meters]
%
%   Outputs:
%   enu = n x 3 matrix of east/north/up data after transformation
%       units are meters
%
%   See also ENU2LLA, LLA2ECR, ECR2ENU, TRANSFORM.

[ecr] = transform('lla2ecr',lla);
[enu] = transform('ecr2enu',ecr,origin_lla);
