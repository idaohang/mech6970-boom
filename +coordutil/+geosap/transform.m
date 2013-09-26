function [output] = transform(type,input,option)
%TRANSFORM Perform coordinate transforms.
%   [output] = transform(type,input,option)
%
%   Inputs:
%   type - string containing the coordinate transform requested:
%       'ecr2lla' 'lla2ecr' 'ecr2enu' 'enu2ecr'
%   input - matrix of data to be transformed
%   option - for enu transformations, contains the desired origin position
%       origin_lla = [latitude longitude altitude] (radians, radians, meters)
%
%   Outputs:
%   output - matrix of transformed data
%
%   See also LLA2ENU.

% note: this error handling is "overbuilt" right now so that it doesn't
% have to be re-tooled if we add support for a bunch more transform types

% make TYPE lowercase so we don't have to ignore case throughout
type = lower(type);

% check that TYPE is supported...
if ~issubfun(type)
    error('GEOSAP:transform','Unsupported transform')
end

% check that input data has 3 columns
[num_rows,num_cols,num_trajs] = size(input);
if ~(num_cols == 3)
    error('GEOSAP:transform','Invalid number of columns')
end

nin = min(num_cols,nargin2(type)-1);
% assume a transform cannot generate more outputs than it has inputs
nout = min(nin,nargout2(type));

% extract the data from the columns
% tried it the num2cell(input,[1 3]) way, but it was slower...
for i = fliplr(1:nin)
    inargs{i} = input(:,i,:);
end

% if option is undefined, make it empty
if nargin < 3
    option = [];
end

% perform the transform
[varargout{1:nout}] = feval(type,option,inargs{:});

% put data back into columns
output = [varargout{:}];



% ------------------------SUBFUNCTIONS--------------------------------
function [lat,lon,h] = ecr2lla(optional,x,y,z)
%ECR2LLA Earth Centered Rotating to Lat/Lon/Alt tranformation.
%   [LAT,LON,H] = ECR2LLA(X,Y,Z);
%
%   See also LLA2ECR.

Re = 6378137;       % Earth_Equatorial_Radius, in m
Rp = 6356752;       % Earth_Polar_Radius, in m
f  = (Re-Rp)/Re;       % Earth_Flattening Coefficient
e2 = sqrt(2*f-f^2)^2;    % Earth_Eccentricity squared
er = e2/(1 - e2);

lon = atan2(y,x);
r = sqrt(x.^2 + y.^2 + z.^2);
sph = z./r;
cph = sqrt(1 - sph.^2);
rt = Re./sqrt(1 + er.*sph.^2);
a = e2.*sph.*cph./(1 - e2.*sph.^2);
b = a.*rt./r;
lat = asin(z./r) + b;
h = (r - rt).*(1 - a.*b./2);



function [x,y,z] = lla2ecr(optional,lat,lon,h)
%LLA2ECR Lat/Lon/Alt to Earth Centered Rotating Transform.
%   [X,Y,Z] = LLA2ECR(LAT,LON,H);
%
%   See also ECR2LLA.

Re = 6378137;       % Earth_Equatorial_Radius, in m
Rp = 6356752;       % Earth_Polar_Radius, in m
f  = (Re-Rp)/Re;    % Earth_Flattening Coefficient

% Compute the ECR coordinates
temp = Re./sqrt(1 + ((1 - f).*tan(lat)).^2);

x = (temp + h.*cos(lat)).*cos(lon);
y = (temp + h.*cos(lat)).*sin(lon);
z = temp.*((1 - f)^2.*tan(lat)) + h.*sin(lat);



function [x,y,z] = ecr2enu(origin_lla,x,y,z)
%ECR2ENU Earth Centered Rotating to East/North/Up transformation.
%
%   See also ENU2ECR.

lat = origin_lla(1);
lon = origin_lla(2);
h   = origin_lla(3);

T = [-sin(lon) cos(lon) 0;
    -sin(lat)*cos(lon) -sin(lat)*sin(lon) cos(lat);
    cos(lat)*cos(lon) cos(lat)*sin(lon) sin(lat)];

radar_ecr = transform('lla2ecr',origin_lla);
radar_rrc = T*radar_ecr';

out = T*[x(:) y(:) z(:)]';
% must reshape each piece of data individually or RESHAPE will mix them up
% create a size matrix for use in reshaping output vectors
% typically, [num_rows 1 num_trajs]
% (this allows us to do a matrix multiply on ND data)
outsize = size(x);
x = reshape(out(1,:),outsize) - radar_rrc(1);
y = reshape(out(2,:),outsize) - radar_rrc(2);
z = reshape(out(3,:),outsize) - radar_rrc(3);

% here's the much slower, but more readable way to do it:
% for i = 1:num_trajs
%     output(:,1:3,i) = (T*input(:,1:3,i)')';
% end
% output(:,1:3,:) = output(:,1:3,:) - radar_rrc;



function [x,y,z] = enu2ecr(origin_lla,x,y,z)
%ENU2ECR East/North/Up to Earth Centered Rotating transformation.
%
%   See also ECR2ENU.

lat = origin_lla(1);
lon = origin_lla(2);
h   = origin_lla(3);

T = [-sin(lon) cos(lon) 0;
    -sin(lat)*cos(lon) -sin(lat)*sin(lon) cos(lat);
    cos(lat)*cos(lon) cos(lat)*sin(lon) sin(lat)];
T = T';

% create a size matrix for use in reshaping output vectors
% typically, [num_rows 1 num_trajs]
% (this allows us to do a matrix multiply on ND data)
outsize = size(x);

radar_ecr = transform('gc2ecr',origin_lla);

out = T*[x(:) y(:) z(:)]';
x = reshape(out(1,:),outsize) + radar_ecr(1);
y = reshape(out(2,:),outsize) + radar_ecr(2);
z = reshape(out(3,:),outsize) + radar_ecr(3);



function [nargs] = nargin2(type)
%NARGIN Overrides the built-in NARGIN, which doesn't work on subfunctions.
%   since NARGIN and NARGOUT don't work on subfunctions, we must do this...
%   (could put the subfunctions in a PRIVATE directory, but NARGIN doesn't work on them, either)
%   (as long as we're doing this, we might as well use VARARGIN and VARARGOUT)
%   (but let's not, for now, just in case NARGIN/NARGOUT get fixed...)
%   (but in Matlab 6.5, it looks like they chose to state that NARGIN only
%   works on separate m-files, rather than fix it for subfunctions...)
%
%   See also NARGIN.

fcnnames = {'ecr2lla' 'lla2ecr' 'ecr2enu' 'enu2ecr'};
fcnin  = [4 4 4 4];
ind = find(strcmp(fcnnames,type));
nargs = fcnin(ind);



function [nargs] = nargout2(type)
%NARGOUT2 Overrides the built-in NARGOUT, which doesn't work on subfunctions.
%   [nargs] = nargout2(type);
%
%   See also NARGOUT.

fcnnames = {'ecr2lla' 'lla2ecr' 'ecr2enu' 'enu2ecr'};
fcnout = [3 3 3 3];
ind = find(strcmp(fcnnames,type));
nargs = fcnout(ind);



function [flag] = issubfun(type)
%ISSUBFUN Since EXIST doesn't tell if something is a subfun, we do this...
%
%   See also ISVARNAME.

fcnnames = {'ecr2lla' 'lla2ecr' 'ecr2enu' 'enu2ecr'};
flag = any(strcmp(fcnnames,type));
