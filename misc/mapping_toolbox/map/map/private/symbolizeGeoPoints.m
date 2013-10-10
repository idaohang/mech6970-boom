function hg = symbolizeGeoPoints(...
                    mstruct, S, symspec, defaultProps, otherProps)
%symbolizeGeoPoints Symbolize and display point geostruct
%
%   hg = symbolizeGeoPoints(mstruct, S, symspec, plotFcn, defaultProps,
%   otherProps) uses the symbol spec in SYMSPEC to map the attributes of
%   the features in the point geostruct S to graphics properties.  Then, it
%   calls the plot function using the function handle given in PLOTFCN.
%   The return value is a handle to an hggroup object containing one
%   graphics object for each feature in S.  Both defaultProps and
%   otherProps are cell arrays containing name-value pairs for graphic
%   properties,
%
%       {prop1, val1, prop2, val2, ...}
%
%   In the event of conflicts, otherProps overrides the symbolization
%   rules from symbolspec, and both override any property values
%   specified in defaultProps.
%
%   Example
%   -------
%   % Create a map of North America.
%   figure
%   worldmap('na');
%
%   % Read the worldcites.
%   cities = shaperead('worldcities', 'UseGeoCoords', true);
%
%   % Create a SymbolSpec to display Chicago and Denver as blue circles.
%   symspec = makesymbolspec('Point', ...
%     {'Name','Chicago','MarkerEdgeColor','blue','Marker','o'}, ...
%     {'Name','Denver','MarkerEdgeColor','blue','Marker','o'});
%
%   % Display all the cities.
%   hg = symbolizeGeoPoints(gcm, cities, symspec, ...
%                          {'MarkerEdgeColor','green'},{});
%
%   See also symbolizeMapVectors.

% Copyright 2006-2008 The MathWorks, Inc.
% $Revision: 1.1.6.2 $  $Date: 2008/06/16 16:47:03 $ 

% Map geostruct attributes to graphics properties.
properties = attributes2properties(symspec, S);

% Trim the point coordinates and project to map x-y.  When a point is
% trimmed away, place a NaN in its place in the vectors x and y, so that
% both of these vectors have the same number of elements as S.
lat = [S.Lat];
lon = [S.Lon];
x = NaN + zeros(size(lat));
y = x;
for k=1:numel(lat)
   [x(k), y(k)] = ...
      feval(mstruct.mapprojection, mstruct, lat(k), lon(k), 'line', 'forward');
end

% Initialize output to empty, in case S is empty.
hg = reshape([],[0 1]);

% Create graphics objects for the remaining features, such that each has
% the hggroup as its parent.  Iterative in reverse so that the order of
% the children of hg matches the order of the features in S,
% establishing a one-to-one correspondence.
for k = numel(S):-1:1
    % Get the graphics properties that are controlled by the symbol spec
    % and combine them with the rest.
    symspecProps = extractprops(properties,k);
    props = [defaultProps symspecProps otherProps];
    
    if isempty(hg)
        % We don't have an hggroup yet
        h = mappoint(x(k), y(k), props{:});
        if ~isempty(h)
            % We've found our first non-empty object. Create an hggroup
            % object with the same parent, then re-parent this object
            % into the group.
            hg = hggroup('Parent', get(h,'Parent'));
            set(h,'Parent',hg)
        end
    else
        % We already have an hggroup, so specify it as the parent
        % of this feature.
        mappoint(x(k), y(k), props{:}, 'Parent', hg);
    end
end
