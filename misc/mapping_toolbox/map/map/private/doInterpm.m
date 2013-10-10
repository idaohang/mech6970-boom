function [lat, lon] = doInterpm(lat,lon,maxdiff,method,units)
% Core computations performed by INTERPM

% Copyright 2006-2007 The MathWorks, Inc.
% $Revision: 1.1.6.3 $  $Date: 2007/06/04 21:12:03 $

%  Ensure column vectors

lat = lat(:);
lon = lon(:);

%  Compute the maximum angular distance between each latitude
%  and longitude pair of points.

dist  = max( [abs(diff(lat))'; abs(diff(lon))'] )';

%  Find angular differences which exceed the maximum allowed

indx  = find( dist > maxdiff);

if ~isempty(indx)
     steps = ceil(dist(indx)/maxdiff);      %  No points added each location
     totalpts = sum(steps)-length(steps);   %  Total points to be added
     lastpt  = length(lat);                 %  Current last point in data set
	 lat(length(lat)+totalpts) = 0;         %  Pre-allocate output memory
	 lon(length(lon)+totalpts) = 0;
end

%  Fill in points where the maximum angular difference is
%  exceeded.  Linearly interpolate points between the identified
%  two end points.

switch method

    case {'gc','rh'}
        %   Interpolate along great circles or rhumb lines. If the
        %   data is crude enough to require interpolation, it's
        %   sufficient to interpolate on a sphere.

        [lat, lon] = toRadians(units, lat, lon);

        for i=length(indx):-1:1
            %  Set the index in the original vectors and compute the
            %  interpolation steps.
            loc = indx(i);
            
            %  interpolation insert
            [latinsert, loninsert] = doTrack(method,...
                lat(loc + [0 1]), lon(loc + [0 1]), ...
                [1 0], steps(i)+1);

            %   strip trailing NaNs inserted by track
            latinsert(isnan(latinsert)) = [];
            loninsert(isnan(loninsert)) = [];

            %   remove starting and ending points to avoid duplication
            %   when inserted
            latinsert = latinsert(2:length(latinsert)-1);
            loninsert = loninsert(2:length(loninsert)-1);

            %  Fill in the interpolated data
            lat=[lat(1:loc); latinsert; lat(loc+1:lastpt)];
            lon=[lon(1:loc); loninsert; lon(loc+1:lastpt)];

            %  Update the last point of the data set.  Note that since
            %  the output memory is pre-allocated, the current last point
            %  of the data set is not equal to the length of the data vector
            lastpt = lastpt + length(latinsert);
        end

        lon = unwrapMultipart(lon);
        [lat, lon] = fromRadians(units, lat, lon);

    otherwise
        % interpolate in a platte carree space

        for i=length(indx):-1:1
            %  Set the index in the original vectors and compute the
            %  interpolation steps.
            loc = indx(i);
            
            %  -1 eliminates double hit at end of interpolation insert
            factors = (1:steps(i)-1)' ; 

            latinsert = ((lat(loc+1)-lat(loc))/steps(i))*factors + lat(loc);
            loninsert = ((lon(loc+1)-lon(loc))/steps(i))*factors + lon(loc);

            %  Fill in the interpolated data
            lat=[lat(1:loc); latinsert; lat(loc+1:lastpt)];
            lon=[lon(1:loc); loninsert; lon(loc+1:lastpt)];

            %  Update the last point of the data set.  Note that since
            %  the output memory is pre-allocated, the current last point
            %  of the data set is not equal to the length of the data vector
            lastpt = lastpt + length(latinsert);
        end

end
