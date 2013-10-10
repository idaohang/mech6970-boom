%GeoRasterReference Reference raster to geographic coordinates
%
%   A GeoRasterReference object encapsulates the relationship between a
%   geographic coordinate system and a system of "intrinsic coordinates"
%   anchored to the columns and rows of a 2-D spatially referenced
%   raster grid or image. The raster must be sampled regularly in
%   latitude and longitude and its columns and rows must be aligned with
%   meridians and parallels, respectively.
%
%   Use the GEORASTERREF function to construct a GeoRasterReference object,
%   or use refvecToGeoRasterReference or refmatToGeoRasterReference to
%   convert an existing referencing vector or matrix.
%
%   GeoRasterReference properties:
%      Latlim - Latitude limits [southern_limit northern_limit]
%      Lonlim - Longitude limits [western_limit eastern_limit]
%      RasterSize - Number of cells or samples in each spatial dimension
%      RasterInterpretation - Controls handling of raster edges
%      AngleUnits - Unit of angle used for angle-valued properties
%      ColumnsStartFrom - Edge where column indexing starts: 'south' or 'north'
%      RowsStartFrom - Edge where row indexing starts: 'west' or 'east'
%
%   GeoRasterReference properties (SetAccess = private):
%      DeltaLat - Change in latitude with respect to intrinsic Y
%      DeltaLon - Change in longitude with respect to intrinsic X
%      RasterExtentInLatitude - Extent in latitude of the full raster
%      RasterExtentInLongitude - Extent in longitude of the full raster
%      XLimIntrinsic - Limits of raster in intrinsic X [xMin xMax]
%      YLimIntrinsic - Limits of raster in intrinsic Y [yMin yMax]
%      CoordinateSystemType - Type of external system (constant: 'geographic')
%
%   GeoRasterReference methods:
%      GeoRasterReference - Construct GeoRasterReference object
%      sizesMatch - True if object and raster or image are size-compatible
%      intrinsicToGeographic - Convert from intrinsic to geographic coordinates
%      geographicToIntrinsic - Convert from geographic to intrinsic coordinates
%      intrinsicYToLatitude  - Convert from intrinsic Y to latitude
%      intrinsicXToLongitude - Convert from intrinsic X to longitude
%      latitudeToIntrinsicY  - Convert from latitude to intrinsic Y
%      longitudeToIntrinsicX - Convert from longitude to intrinsic X
%      geographicToSub - Geographic coordinates to row and column subscripts
%      contains - True if raster contains latitude-longitude points
%      worldFileMatrix - World file parameters for transformation
%
%   See also GEORASTERREF, spatialref.MapRasterReference, refmatToGeoRasterReference, refvecToGeoRasterReference

% Copyright 2010-2011 The MathWorks, Inc.
% $Revision: 1.1.6.4.2.2 $  $Date: 2011/02/05 19:22:33 $

% Hidden properties (with SetAccess = private)
%
%    FirstCornerLat - Latitude of the (1,1) corner of the raster
%    FirstCornerLon - Longitude of the (1,1) corner of the raster
%    DeltaLatNumerator   - Numerator of rational DeltaLat property
%    DeltaLatDenominator - Denominator of rational DeltaLat property
%    DeltaLonNumerator   - Numerator of rational DeltaLon property
%    DeltaLonDenominator - Denominator of rational DeltaLon property

classdef (Sealed = true) GeoRasterReference
    
    %------------------- Properties: Public + visible --------------------
    
    properties (Dependent = true)
        %Latlim  Latitude limits
        %
        %   Latlim specifies the limits in latitude of the geographic
        %   quadrangle bounding the georeferenced raster.  It is a
        %   two-element vector of the form:
        %
        %             [southern_limit northern_limit]
        Latlim
        
        %Lonlim  Longitude limits
        %
        %   Lonlim specifies the limits in longitude of the geographic
        %   quadrangle bounding the georeferenced raster.  It is a
        %   two-element vector of the form:
        %
        %             [western_limit eastern_limit]
        Lonlim
        
        %RasterSize Number of cells or samples in each spatial dimension
        %
        %   RasterSize is a two-element vector [M N] specifying the
        %   number of rows (M) and columns (N) of the raster or image
        %   associated with the referencing object. In addition, for
        %   convenience, you may assign a size vector having more than
        %   two elements to RasterSize. This flexibility enables
        %   assignments like R.RasterSize = size(RGB), for example, where
        %   RGB is M-by-N-by-3. However, in such cases, only the first two
        %   elements of the size vector will actually be stored. The higher
        %   (non-spatial) dimensions will be ignored. M and N must be
        %   positive in all cases and must be 2 or greater when
        %   RasterInterpretation is 'postings'.
        RasterSize
        
        %RasterInterpretation Controls handling of raster edges
        %
        %   RasterInterpretation is a string that equals 'cells' or
        %   'postings'.
        RasterInterpretation
        
        %AngleUnits Unit of angle used for angle-valued properties
        %
        %   AngleUnits is a string that equals 'degrees'.
        AngleUnits
        
        %ColumnsStartFrom Edge from which column indexing starts
        %
        %   ColumnsStartFrom is a string that equals 'south' or 'north'.
        ColumnsStartFrom
        
        %RowsStartFrom Edge from which row indexing starts
        %
        %   RowsStartFrom is a string that equals 'west' or 'east'.
        RowsStartFrom
    end
    
    properties (Dependent = true, SetAccess = private)
        % DeltaLat - Change in latitude with respect to intrinsic Y
        %
        %    DeltaLat is the amount by which latitude increases or
        %    decreases with respect to an increase of one unit in
        %    intrinsic Y. It is positive when columns start from south,
        %    and negative when columns start from north. Its absolute
        %    value equals the latitude extent of a single cell (when
        %    RasterInterpretation is 'cells') or the latitude separation
        %    of adjacent sample points (when the RasterInterpretation is
        %    'postings').
        DeltaLat
        
        % DeltaLon - Change in longitude with respect to intrinsic X
        %
        %    DeltaLon is the amount by which longitude increases or
        %    decreases with respect to an increase of one unit in
        %    intrinsic X. It is positive when rows start from west and
        %    negative when rows start from east. Its absolute value
        %    equals the longitude extent of a single cell (when
        %    RasterInterpretation is 'cells') or the longitude
        %    separation of adjacent sample points (when the
        %    RasterInterpretation is 'postings').
        DeltaLon
        
        % RasterExtentInLatitude - Extent in latitude of the full raster
        %
        %    RasterExtentInLatitude is the latitude extent ("height") of the
        %    quadrangle covered by the raster.
        RasterExtentInLatitude
        
        % RasterExtentInLongitude - Extent in longitude of the full raster
        %
        %    RasterExtentInLongitude is the longitude extent ("width") of the
        %    quadrangle covered by the raster.
        RasterExtentInLongitude
    end
    
    properties (SetAccess = private, Transient = true)
        % XLimIntrinsic - Limits of raster in intrinsic X [xMin xMax]
        %
        %    XLimIntrinsic is a two-element row vector. For an M-by-N
        %    raster with RasterInterpretation equal to 'postings' it
        %    equals [1 N], and for 'cells' it equals [0.5, N + 0.5].
        XLimIntrinsic;
        
        % YLimIntrinsic - Limits of raster in intrinsic Y [yMin yMax]
        %
        %    YLimIntrinsic is a two-element row vector. For an M-by-N
        %    raster with RasterInterpretation equal to 'postings' it
        %    equals [1 M], and for 'cells' it equals [0.5, M + 0.5].
        YLimIntrinsic;
    end
    
    properties (Constant = true)
        % CoordinateSystemType - Type of external system (constant: 'geographic')
        %
        %   CoordinateSystemType describes the type of coordinate system
        %   represented to which the image or raster is referenced. It
        %   is a constant string with value 'geographic'.
        CoordinateSystemType = 'geographic';
    end
    
    properties (SetAccess = private, Hidden = true)
        %FirstCornerLat - Latitude of the (1,1) corner of the raster
        %
        %   Latitude of the outermost corner of the first cell (1,1) of
        %   the raster (if RasterInterpretation is 'cells') or the first
        %   sample point (if RasterInterpretation is 'postings').
        FirstCornerLat = 0.5;
        
        %FirstCornerLon - Longitude of the (1,1) corner of the raster
        %
        %   Longitude of the outermost corner of the first cell (1,1) of
        %   the raster (if RasterInterpretation is 'cells') or the first
        %   sample point (if RasterInterpretation is 'postings').
        FirstCornerLon = 0.5;
        
        %DeltaLatNumerator - Numerator of rational DeltaLat property
        %
        %   DeltaLatNumerator is a non-negative real number which, when
        %   divided by DeltaLatDenominator, defines a signed north-south
        %   cell size (when RasterInterpretation is 'cells') or sample
        %   spacing (when RasterInterpretation is 'postings').  A
        %   positive value indicates that columns run from south to
        %   north, whereas a negative value indicates that columns run
        %   from north to south.
        DeltaLatNumerator = 1;
        
        %DeltaLatDenominator - Denominator of rational DeltaLat property
        %
        %   DeltaLatDenominator is a (strictly) positive real number
        %   which, when divided into DeltaLatNumerator, defines a signed
        %   north-south cell size (when RasterInterpretation is 'cells')
        %   or sample spacing (when RasterInterpretation is 'postings').
        DeltaLatDenominator = 1;
        
        %DeltaLonNumerator - Numerator of rational DeltaLon property
        %
        %   DeltaLonNumerator is a non-negative real number which, when
        %   divided by DeltaLonDenominator, defines a signed east-west
        %   cell size (when RasterInterpretation is 'cells') or sample
        %   spacing (when RasterInterpretation is 'postings').  A
        %   positive value indicates that rows run from west to
        %   east, whereas a negative value indicates that rows run
        %   from east to west.
        DeltaLonNumerator = 1;
        
        %DeltaLonDenominator - Denominator rational DeltaLon property
        %
        %   DeltaLonDenominator is a (strictly) positive real number
        %   which, when divided into DeltaLonNumerator, defines a signed
        %   east-west cell size (when RasterInterpretation is 'cells')
        %   or sample spacing (when RasterInterpretation is 'postings').
        DeltaLonDenominator = 1;
    end
    
    %---------------- Properties: Private + hidden ---------------------
    
    properties (Access = private, Hidden = true)
        % The RasterSize and RasterInterpretation properties depend solely
        % on a hidden instance of an IntrinsicRaster2D object. It supports
        % the public sizesMatch(), geographicToSub(), and contains()
        % methods, as well as the get methods for Latlim and Lonlim and
        % several private methods.
        Intrinsic = spatialref.IntrinsicRaster2D();
        
        % A place to hold the value of the public AngleUnits property.
        pAngleUnits = 'degrees';
    end
    
    %-------------------------- load object ----------------------------
    
    methods (Static)
        
        function self = loadobj(self)
            % These two transient properties are not explicitly dependent,
            % but they are architecturally dependent, and need to be set.
            self.XLimIntrinsic = self.Intrinsic.XLim;
            self.YLimIntrinsic = self.Intrinsic.YLim;
        end
        
    end
    
    %-------------- Constructor and ordinary methods -------------------
    
    methods
        
        function self = GeoRasterReference(varargin)
            %GeoRasterReference Construct GeoRasterReference object
            %
            %   R = spatialref.GeoRasterReference() constructs a
            %   GeoRasterReference object with the following default
            %   property settings:
            %
            %                         Latlim: [0.5 2.5]
            %                         Lonlim: [0.5 2.5]
            %                     RasterSize: [2 2]
            %           RasterInterpretation: 'cells'
            %                     AngleUnits: 'degrees'
            %               ColumnsStartFrom: 'south'
            %                  RowsStartFrom: 'west'
            %                       DeltaLat: 1
            %                       DeltaLon: 1
            %         RasterExtentInLatitude: 2
            %        RasterExtentInLongitude: 2
            %                  XLimIntrinsic: [0.5 2.5]
            %                  YLimIntrinsic: [0.5 2.5]
            %           CoordinateSystemType: 'geographic'
            %
            %   R = spatialref.GeoRasterReference( ...
            %       rasterSize, rasterInterpretation, angleUnits, ...
            %       firstCornerLat, firstCornerLon, ...
            %       deltaLatNumerator, deltaLatDenominator, ...
            %       deltaLonNumerator, deltaLonDenominator)
            %   is a special advanced syntax for constructing a
            %   GeoRasterReference object from the following inputs (all
            %   nine must be provided when the constructor is called):
            %
            %     rasterSize -- A valid MATLAB size vector (as returned
            %        by SIZE) corresponding to the size of a raster or
            %        image to be used in conjunction with the
            %        IntrinsicRaster2D object. rasterSize must have at
            %        least two elements, but may have more. In this case
            %        only the first two (corresponding to the "spatial
            %        dimensions") are used and the others are ignored.
            %        For example, if rasterSize = size(RGB) where RBG is
            %        an RGB image, the rasterSize will have the form
            %        [M N 3] and the RasterSize property will be set to
            %        [M N].  M and N must be positive in all cases and
            %        must be 2 or greater when RasterInterpretation is
            %        'postings'.
            %
            %     rasterInterpretation -- A string: either 'cells' or
            %     'postings'
            %
            %     angleUnits -- A string (must equal 'degrees')
            %
            %     firstCornerLat, firstCornerLon -- Scalar values
            %        defining the latitude and longitude position of the
            %        outermost corner of the first cell (1,1) of the
            %        raster (if rasterInterpretation is 'cells') or the
            %        first sample point (if rasterInterpretation is
            %        'postings')
            %
            %     deltaLatNumerator   -- Nonzero real number
            %     deltaLatDenominator -- Positive real number
            %
            %        The ratio deltaLatNumerator/deltaLatDenominator
            %        defines a signed north-south cell size (when
            %        rasterInterpretation is 'cells') or sample spacing
            %        (when rasterInterpretation is 'postings').  A
            %        positive value indicates that columns run from
            %        south to north, whereas a negative value indicates
            %        that columns run from north to south.
            %
            %     deltaLonNumerator   -- Nonzero real number
            %     deltaLonDenominator -- Positive real number
            %
            %        The ratio deltaLonNumerator/deltaLonDenominator
            %        defines a signed east-west cell size (when
            %        rasterInterpretation is 'cells') or sample spacing
            %        (when rasterInterpretation is 'postings').  A
            %        positive value indicates that rows run from west to
            %        east, whereas a negative value indicates that rows
            %        run from east to west.
            %
            %   Examples using advanced syntax
            %   ------------------------------
            %   % Construct a referencing object for a global raster
            %   % comprising 180-by-360 one-degree cells, with rows that
            %   % start at longitude -180, and with the first cell
            %   % located in the northwest corner.
            %   R = spatialref.GeoRasterReference( ...
            %       [180 360], 'cells', 'degrees', 90, -180, -1, 1, 1, 1)
            %
            %   % Construct a referencing object for the DTED Level 0
            %   % file that includes Sagarmatha (Mount Everest).
            %   R = spatialref.GeoRasterReference([121 121], ...
            %       'postings', 'degrees', 27, 86, 1, 120, 1, 120)
            %
            %   See also GEORASTERREF.
            
            if nargin == 0
                % Default constructor.
                
                % These two properties are not explicitly dependent, but
                % they are architecturally dependent, and don't have
                % default values -- so assign them now.
                self.XLimIntrinsic = self.Intrinsic.XLim;
                self.YLimIntrinsic = self.Intrinsic.YLim;
            elseif nargin == 9
                % Construct from complete list of defining properties.
                
                [rasterSize, rasterInterpretation, angleUnits, ...
                    firstCornerLat, firstCornerLon, ...
                    deltaLatNumerator, deltaLatDenominator, ...
                    deltaLonNumerator, deltaLonDenominator] = deal(varargin{:});
                
                % Construct (and validate) hidden component object.
                self.Intrinsic = spatialref.IntrinsicRaster2D( ...
                    rasterSize, rasterInterpretation);                
                self.Intrinsic.validate()
                
                % Note: It's essential to assign AngleUnits before
                % setting any of the angle-valued properties.
                self.AngleUnits = angleUnits;
                
                % Make sure to set XLimIntrinsic and YLimIntrinsic
                % before setting properties related to latitude and
                % longitude, so they can be used in their validation.
                self.XLimIntrinsic = self.Intrinsic.XLim;
                self.YLimIntrinsic = self.Intrinsic.YLim;
                
                self = self.setLatitudeProperties( ...
                    firstCornerLat, deltaLatNumerator, deltaLatDenominator);
                
                self = self.setLongitudeProperties( ...
                    firstCornerLon, deltaLonNumerator, deltaLonDenominator);
            else
                error('map:spatialref:expectedNineInputs', ...
                    'Either 0 or exactly 9 input arguments are required.')
            end
        end
        
        
        function tf = sizesMatch(self,A)
            %sizesMatch True if object and raster or image are size-compatible
            %
            %   TF = R.sizesMatch(A) returns true if the size of the raster
            %   (or image) A is consistent with the RasterSize property of
            %   the referencing object R. That is,
            %
            %           R.RasterSize == [size(A,1) size(A,2)].
            
            tf = self.Intrinsic.sizesMatch(A);
        end
        
        
        function [lat, lon] = intrinsicToGeographic(self, xi, yi)
            %intrinsicToGeographic Convert from intrinsic to geographic coordinates
            %
            %   [LAT, LON] = R.intrinsicToGeographic(xIntrinsic, yIntrinsic)
            %   returns the geographic coordinates (LAT, LON) of a set
            %   of points given their intrinsic coordinates (xIntrinsic,
            %   yIntrinsic), based on the relationship defined by the
            %   referencing object R. xIntrinsic and yIntrinsic must
            %   have the same size. LAT and LON will have the same size
            %   as xIntrinsic and yIntrinsic. The input may include
            %   points that fall outside the limits of the raster (or
            %   image). Latitudes and longitudes for such points are
            %   linearly extrapolated outside the geographic quadrangle
            %   bounding the raster, but for any point that extrapolates
            %   to a latitude beyond the poles (latitude < -90 degrees
            %   or latitude > 90 degrees), the values of both LAT and
            %   LON are set to NaN.
            
            validateCoordinatePairs(xi, yi, ...
                'GeoRasterReference.intrinsicToGeographic', ...
                'IntrinsicX', 'IntrinsicY')
            
            % Apply linear scale and shift
            lat = self.intrinsicYToLatitude(yi);
            lon = self.intrinsicXToLongitude(xi);
            
            % Extrapolation beyond the poles causes NaNs to be placed
            % in LAT. Replicate such NaNs in LON also.
            lon(isnan(lat)) = NaN;
        end
        
        
        function [xi, yi] = geographicToIntrinsic(self, lat, lon)
            %geographicToIntrinsic Convert from geographic to intrinsic coordinates
            %
            %   [xIntrinsic, yIntrinsic] = R.geographicToIntrinsic(LAT, LON)
            %   returns the intrinsic coordinates (xIntrinsic, yIntrinsic)
            %   of a set of points given their geographic coordinates
            %   (LAT, LON), based on the relationship defined by the
            %   referencing object R. LAT and LON must have the same
            %   size, and all (non-NaN) elements of LAT must fall within
            %   the interval [-90 90] degrees. xIntrinsic and yIntrinsic
            %   will have the same size as LAT and LON. The input may
            %   include points that fall outside the geographic
            %   quadrangle bounding the raster. As long as their
            %   latitudes are valid, the locations of such points will
            %   be extrapolated outside the bounds of the raster in the
            %   intrinsic coordinate system.
            
            validateCoordinatePairs(lat, lon, ...
                'GeoRasterReference.geographicToIntrinsic', 'LAT', 'LON')
            
            xi = self.longitudeToIntrinsicX(lon);
            yi = self.latitudeToIntrinsicY(lat);
            
            % Latitudes beyond the poles cause NaNs to be placed
            % in yi. Replicate such NaNs in xi also.
            xi(isnan(yi)) = NaN;
        end
        
        
        function lat = intrinsicYToLatitude(self, yi)
            %intrinsicYToLatitude Convert from intrinsic Y to latitude
            %
            %   LAT = R.intrinsicYToLatitude(yIntrinsic) returns the
            %   latitude of the small circle corresponding to the line
            %   y = yIntrinsic, based on the relationship defined by the
            %   referencing object R. The input may include values that
            %   fall completely outside the intrinsic Y-limits of the
            %   raster (or image). In this case latitude is extrapolated
            %   outside the latitude limits, but for input values that
            %   extrapolate to latitudes beyond the poles (latitude < -90
            %   degrees or latitude > 90 degrees), the value of LAT is
            %   set to NaN. NaN-valued elements of yIntrinsic map to
            %   NaNs in LAT.
            
            lat = self.FirstCornerLat + (yi - self.YLimIntrinsic(1)) ...
                .* self.DeltaLatNumerator ./ self.DeltaLatDenominator;
            
            lat(beyondPole(lat, self.pAngleUnits)) = NaN;
        end
        
        
        function lon = intrinsicXToLongitude(self, xi)
            %intrinsicXToLongitude Convert from intrinsic X to longitude
            %
            %   LON = R.intrinsicXToLongitude(xIntrinsic) returns the
            %   longitude of the meridian corresponding to the line
            %   x = xIntrinsic, based on the relationship defined by the
            %   referencing object R. The input may include values that
            %   fall completely outside the intrinsic X-limits of the
            %   raster (or image). In this case, longitude is
            %   extrapolated outside the longitude limits. NaN-valued
            %   elements of xIntrinsic map to NaNs in LON.
            
            lon = self.FirstCornerLon + (xi - self.XLimIntrinsic(1)) ...
                .* self.DeltaLonNumerator ./ self.DeltaLonDenominator;
        end
        
        
        function yi = latitudeToIntrinsicY(self, lat)
            %latitudeToIntrinsicY Convert from latitude to intrinsic Y
            %
            %   yIntrinsic = R.latitudeToIntrinsicY(LAT) returns the
            %   intrinsic Y value of the line corresponding to the small
            %   circle at latitude LAT, based on the relationship
            %   defined by the referencing object R. The input may
            %   include values that fall completely outside the latitude
            %   limits of the raster (or image). In this case yIntrinsic
            %   is either extrapolated outside the intrinsic Y limits --
            %   for elements of LAT that fall within the interval
            %   [-90 90] degrees, or set to NaN -- for elements of LAT
            %   that do not correspond to valid latitudes. NaN-valued
            %   elements of LAT map to NaNs in yIntrinsic.
            
            % Elements of LAT are less than -90 degrees or
            % that exceed +90 degrees should map to NaN.
            lat(beyondPole(lat, self.pAngleUnits)) = NaN;
            
            % Shift and scale latitude
            yi = self.YLimIntrinsic(1) + (lat - self.FirstCornerLat) ...
                .* self.DeltaLatDenominator ./ self.DeltaLatNumerator;
        end
        
        
        function xi = longitudeToIntrinsicX(self, lon)
            %longitudeToIntrinsicX Convert from longitude to intrinsic X
            %
            %   xIntrinsic = R.longitudeToIntrinsicX(LON) returns the
            %   intrinsic X value of the line corresponding to the
            %   meridian at longitude LON, based on the relationship
            %   defined by the referencing object R. The input may
            %   include values that fall completely outside the
            %   longitude limits of the raster (or image). In this case
            %   xIntrinsic is extrapolated outside the intrinsic X
            %   limits. NaN-valued elements of LON map to NaNs in
            %   xIntrinsic.
            
            lonlim = self.getLonlim();
            w = lonlim(1);
            e = lonlim(2);
            
            % Adjust longitude wrapping to get within the limits,
            % whenever possible.
            angleUnits = self.pAngleUnits;
            if (e - w) <= fullCycle(angleUnits)
                rowsRunWestToEast = (self.DeltaLonNumerator > 0);
                if rowsRunWestToEast
                    % Wrap to interval self.FirstCornerLon + [0 360]
                    lon = w + wrapToCycle(lon - w, angleUnits);
                else
                    % Wrap to interval self.FirstCornerLon + [-360 0]
                    lon = e - wrapToCycle(e - lon, angleUnits);
                end
            else
                % Any longitude can be wrapped to fall within the
                % interval [w e], and in fact there's more than one
                % solution for certain longitudes. Resolve the ambiguity
                % by moving longitudes that are west of the western
                % limit the minimal number of cycles to the east that
                % puts them within the limits. Likewise, move longitudes
                % that exceed the eastern limit the minimum number of
                % cycles to the west.
                offToWest = lon < w;
                lon(offToWest) = ...
                    w + wrapToCycle(lon(offToWest) - w, angleUnits);
                
                offToEast = lon > e;
                t = e - fullCycle(angleUnits);
                lon(offToEast) ...
                    = t + wrapToCycle(lon(offToEast) - t, angleUnits);
            end
            
            % Shift and scale longitude
            xi = self.XLimIntrinsic(1) + (lon - self.FirstCornerLon) ...
                .* self.DeltaLonDenominator ./ self.DeltaLonNumerator;
        end
        
        
        function [row,col] = geographicToSub(self, lat, lon)
            %geographicToSub Geographic coordinates to row and column subscripts
            %
            %   [I,J] = R.geographicToSub(LAT,LON) returns the subscript
            %   arrays I and J. When the referencing object R has
            %   RasterInterpretation 'cells', these are the row and column
            %   subscripts of the raster cells (or image pixels) containing
            %   each element of a set of points given their geographic
            %   coordinates (LAT, LON). If R.RasterInterpretation is
            %   'postings', then the subscripts refer to the nearest sample
            %   point (posting). LAT and LON must have the same size. I and
            %   J will have the same size as LAT and LON. For an M-by-N
            %   raster, 1 <= I <= M and 1 <= J <= N, except when a point
            %   LAT(k),LON(k) falls outside the image, as defined by
            %   R.contains(lat, lon), then both I(k) and J(k) are NaN.
            
            % Note: geographicToIntrinsic validates lat and lon
            [xi, yi] = self.geographicToIntrinsic(lat, lon);
            [row, col] = self.Intrinsic.intrinsicToSub(xi, yi);
        end
        
        
        function tf = contains(self, lat, lon)
            %contains True if raster contains latitude-longitude points
            %
            %   TF = R.contains(LAT,LON) returns a logical array TF
            %   having the same size as LAT and LON such that TF(k) is
            %   true if and only if the point (LAT(k),LON(k)) falls
            %   within the bounds of the raster associated with
            %   referencing object R. Elements of LON can be wrapped
            %   arbitrarily without affecting the result.
            
            % Note: This implementation is a minor adaptation of the
            % Mapping Toolbox function INGEOQUAD, in which we simply
            % omit the wrapping step when computing "londiff".
            % And it's generalized to work with radians as well as
            % degrees.
            
            validateCoordinatePairs(lat, lon, ...
                'GeoRasterReference.geographicToIntrinsic', 'LAT', 'LON')
            
            latlim = self.getLatlim();
            lonlim = self.getLonlim();
            
            % Initialize to include all points.
            tf = true(size(lat));
            
            % Eliminate points that fall outside the latitude limits.
            inlatzone = (latlim(1) <= lat) & (lat <= latlim(2));
            tf(~inlatzone) = false;
            
            % Eliminate points that fall outside the longitude limits.
            londiff = lonlim(2) - lonlim(1);  % No need to wrap here
            inlonzone = (wrapToCycle( ...
                lon - lonlim(1), self.pAngleUnits) <= londiff);
            tf(~inlonzone) = false;
        end
        
        
        function W = worldFileMatrix(self)
            %worldFileMatrix - World file parameters for transformation
            %
            %   W = R.worldFileMatrix returns a 2-by-3 world file matrix.
            %   Each of the 6 elements in W matches one of the lines in a
            %   world file corresponding to the transformation defined by
            %   the referencing object R.
            %
            %   Given W with the form:
            %
            %                    W = [A B C;
            %                         D E F],
            %
            %   a point (xi, yi) in intrinsic coordinates maps to a point
            %   (lat, lon) in geographic coordinates like this:
            %
            %         lon = A * (xi - 1) + B * (yi - 1) + C
            %         lat = D * (xi - 1) + E * (yi - 1) + F
            %
            %   or, more compactly, [lon lat]' = W * [(xi - 1) (yi - 1) 1]'.
            %   The -1s allow the world file matrix to work with the
            %   Mapping Toolbox convention for intrinsic coordinates, which
            %   is consistent with the 1-based indexing used throughout
            %   MATLAB. W is stored in a world file with one term per line
            %   in column-major order: A, D, B, E, C, F.  That is, a world
            %   file contains the elements of W in the following order:
            % 
            %         W(1,1)
            %         W(2,1)
            %         W(1,2)
            %         W(2,2)
            %         W(1,3)
            %         W(2,3).
            %
            %   The expressions above hold for a general affine
            %   transformation, but in the matrix returned by this method
            %   B, D, W(2,1), and W(1,2) are identically 0 because
            %   longitude depends only on intrinsic X and latitude depends
            %   only on intrinsic Y.
            %
            %   See also WORLDFILEREAD, WORLDFILEWRITE.
            
            if strcmp(self.RasterInterpretation,'postings')
                c = self.FirstCornerLon;
                f = self.FirstCornerLat;
            else
                c = self.intrinsicXToLongitude(1);
                f = self.intrinsicYToLatitude(1);
            end
            W = [self.DeltaLon         0          c;
                       0        self.DeltaLat     f];
        end
        
    end
    
    %------------------------- Overloaded disp --------------------------
    
    methods
        
        function disp(self)
            % Override the default to display the values of the DeltaLat
            % and DeltaLon properties as rational numbers.
            
            s = evalc('builtin(''disp'',self)');
            
            if isscalar(self)
                if self.DeltaLatDenominator ~= 1
                    s = replaceValueWithRatio(s, 'DeltaLat', ...
                        self.DeltaLatNumerator, ...
                        self.DeltaLatDenominator);
                end
                
                if self.DeltaLonDenominator ~= 1
                    s = replaceValueWithRatio(s, 'DeltaLon', ...
                        self.DeltaLonNumerator, ...
                        self.DeltaLonDenominator);
                end
            end
            
            fprintf('%s',s)
        end
        
    end
    
    %-------------------------- Set methods ----------------------------
    
    methods
        
        function self = set.RasterSize(self, rasterSize)           
            
            % Save current values of dependent properties.
            latlim = self.getLatlim();
            lonlim = self.getLonlim();
            
            % Update (and validate) RasterSize property.
            self.Intrinsic.RasterSize = rasterSize;
            self.Intrinsic.validate()
            
            % Update cached values (right away).
            self.XLimIntrinsic = self.Intrinsic.XLim;
            self.YLimIntrinsic = self.Intrinsic.YLim;
            
            % Reset "delta" properties while maintaining latitude and
            % longitude limits, and column and row directions.
            self = constrainToFitLatlim(self, latlim);
            self = constrainToFitLonlim(self, lonlim);
            
        end
        
        
        function self = set.RasterInterpretation(self, rasterInterpretation)
            
            % Save current values of dependent properties.
            latlim = self.getLatlim();
            lonlim = self.getLonlim();
            
            % Update (and validate) RasterInterpretation property.
            self.Intrinsic.RasterInterpretation = rasterInterpretation;
            self.Intrinsic.validate()
            
            % Update cached values (right away).
            self.XLimIntrinsic = self.Intrinsic.XLim;
            self.YLimIntrinsic = self.Intrinsic.YLim;
            
            % Reset "delta" properties while maintaining latitude and
            % longitude limits, and column and row directions.
            self = constrainToFitLatlim(self, latlim);
            self = constrainToFitLonlim(self, lonlim);
        end
        
        
        function self = set.AngleUnits(self, angleUnits)
            angleUnits = validatestring(angleUnits, {'degrees'});
            self.pAngleUnits = angleUnits;
        end
        
        
        function self = set.Latlim(self, latlim)
            validateattributes(latlim, ...
                {'double'}, {'real','row','finite','size',[1 2]}, ...
                'spatialref.GeoRasterReference.set.Latlim', 'latlim')
            
            assert(~any(beyondPole(latlim, self.pAngleUnits)), ...
                'map:spatialref:invalidLatlim', ...
                'Latitude limits must be within [-90 90] degrees.')
            
            assert(latlim(1) < latlim(2), ...
                'map:spatialref:expectedAscendingLimits', ...
                'Elements of %s must be ascending in value.', ...
                'latlim')
            
            % Reset delta latitude properties while maintaining raster
            % size and column direction.
            self = constrainToFitLatlim(self, latlim);            
        end
        
        
        function self = set.Lonlim(self, lonlim)
            validateattributes(lonlim, ...
                {'double'}, {'real','row','finite','size',[1 2]}, ...
                'spatialref.GeoRasterReference.set.Lonlim', 'lonlim')
            
            assert(lonlim(1) < lonlim(2), ...
                'map:spatialref:expectedAscendingLimits', ...
                'Elements of %s must be ascending in value.', ...
                'lonlim')
            
            % Reset delta longitude properties while maintaining raster
            % size and row direction.
            self = constrainToFitLonlim(self, lonlim);
        end
        
        
        function self = set.ColumnsStartFrom(self, edge)
            edge = validatestring(edge,{'north','south'});
            latlim = self.getLatlim();
            if strcmp(edge,'south')
                % Columns run south to north
                self.FirstCornerLat = latlim(1);
                self.DeltaLatNumerator = abs(self.DeltaLatNumerator);
            else
                % Columns run north to south
                self.FirstCornerLat = latlim(2);
                self.DeltaLatNumerator = -abs(self.DeltaLatNumerator);
            end
        end
        
        
        function self = set.RowsStartFrom(self, edge)
            edge = validatestring(edge,{'east','west'});
            lonlim = self.getLonlim();
            if strcmp(edge,'west')
                % Rows run west to east
                self.FirstCornerLon = lonlim(1);
                self.DeltaLonNumerator = abs(self.DeltaLonNumerator);
            else
                % Rows run east to west
                self.FirstCornerLon = lonlim(2);
                self.DeltaLonNumerator = -abs(self.DeltaLonNumerator);
            end
        end
        
    end
    
    %----------------- Get methods for public properties ------------------
    
    methods
        
        function rasterSize = get.RasterSize(self)
            rasterSize = self.Intrinsic.RasterSize;
        end
        
        
        function rasterInterp = get.RasterInterpretation(self)
            rasterInterp = self.Intrinsic.RasterInterpretation;
        end
        
        
        function angleUnits = get.AngleUnits(self)
            angleUnits = self.pAngleUnits;
        end
        
        
        function limits = get.Latlim(self)
            limits = self.getLatlim();
        end
        
        
        function limits = get.Lonlim(self)
            limits = self.getLonlim();
        end
        
        
        function edge = get.ColumnsStartFrom(self)
            if self.DeltaLatNumerator > 0
                edge = 'south';
            else
                edge = 'north';
            end
        end
        
        
        function edge = get.RowsStartFrom(self)
            if self.DeltaLonNumerator > 0
                edge = 'west';
            else
                edge = 'east';
            end
        end
        
        
        function delta = get.DeltaLat(self)
            delta = self.DeltaLatNumerator / self.DeltaLatDenominator;
        end
        
        
        function delta = get.DeltaLon(self)
            delta = self.DeltaLonNumerator / self.DeltaLonDenominator;
        end
        
        
        function extent = get.RasterExtentInLatitude(self)
            extent = diff(self.getLatlim);
        end
        
        
        function extent = get.RasterExtentInLongitude(self)
            rasterExtentInXIntrinsic = abs(diff(self.XLimIntrinsic));
            extent = abs(self.DeltaLonNumerator) ...
                * rasterExtentInXIntrinsic / self.DeltaLonDenominator;
        end
                
    end
    
    %---------------------------- Private methods ------------------------
    
    methods (Access = private)
        
        function latlim = getLatlim(self)
            yi = self.YLimIntrinsic;
            columnsRunSouthToNorth = (self.DeltaLatNumerator > 0);
            if columnsRunSouthToNorth
                latlim = self.intrinsicYToLatitude(yi);
            else
                latlim = self.intrinsicYToLatitude(yi([2 1]));
            end
        end
        
        
        function lonlim = getLonlim(self)
            xi = self.XLimIntrinsic;
            rowsRunWestToEast = (self.DeltaLonNumerator > 0);
            if rowsRunWestToEast
                lonlim = self.intrinsicXToLongitude(xi);
            else
                lonlim = self.intrinsicXToLongitude(xi([2 1]));
            end
        end
        
        
        function self = constrainToFitLatlim(self, latlim)
            % Constrain properties controlling latitude vs. intrinsicY
            %
            % If the raster has one or more rows, set the values of
            % these defining properties:
            %
            %     FirstCornerLat
            %     DeltaLatNumerator
            %     DeltaLatDenominator
            %
            % to be consistent with the specified LATLIM value, unless
            % the RasterInterpretation is 'postings' and the raster
            % has only one row. The extent in latitude (and intrinsic Y) of
            % such a raster is 0, so it requires special handling.
            
            s = latlim(1);
            n = latlim(2);
            if self.DeltaLatNumerator > 0
                dlat = n - s;
                self.FirstCornerLat = s;
            else
                dlat = s - n;
                self.FirstCornerLat = n;
            end
            [self.DeltaLatNumerator, self.DeltaLatDenominator] ...
                = simplifyRatio(dlat, diff(self.YLimIntrinsic));
        end
        
        
        function self = setLatitudeProperties(self, ...
                firstCornerLat, deltaLatNumerator, deltaLatDenominator)
            % Set properties controlling latitude vs. intrinsicY
            %
            % Set the following properties as a group:
            %
            %      FirstCornerLat
            %      DeltaLatNumerator
            %      DeltaLatDenominator
            %
            % These properties in combination with RasterSize(1)
            % determine the latitude limits and must be validated
            % together.
            
            % Validate individual inputs
            fname = 'GeoRasterReference.setLatitudeProperties';
            
            validateattributes(firstCornerLat, {'double'}, ...
                {'real','scalar','finite'}, fname, 'firstCornerLat')
            
            validateattributes(deltaLatNumerator, {'double'}, ...
                {'real','scalar','finite','nonzero'}, ...
                fname, 'deltaLatNumerator')
            
            validateattributes(deltaLatDenominator, {'double'}, ...
                {'real','scalar','finite','positive'}, ...
                fname, 'deltaLatDenominator')
            
            % Assign property values
            self.FirstCornerLat = firstCornerLat;
            [self.DeltaLatNumerator, self.DeltaLatDenominator] ...
                = simplifyRatio(deltaLatNumerator, deltaLatDenominator);
            
            % Note: At this point the value object SELF has been
            % modified and could have invalid latitude limits, but
            % that's OK because we're about to check them. If there's
            % a problem this method will end in an error rather than
            % returning the modified version of SELF.
            
            % Determine the latitude limits implied by the inputs
            latlim = self.getLatlim();
            
            % Validate the implied latitude limits
            assert(all(isfinite(latlim)), ...
                'map:spatialref:invalidLatProps', ...
                ['In combination with specified number of rows, %d,', ...
                ' the values specified for the %s, %s, and %s properties' ...
                ' imply latitude limits that extend outside the interval' ...
                ' [-90 90] degrees.'], self.RasterSize(1), ...
                'FirstCornerLat', 'DeltaLatNumerator', 'DeltaLatDenominator')
        end
        
        
        function self = constrainToFitLonlim(self, lonlim)
            % Constrain properties controlling longitude vs. intrinsicX
            %
            % If the raster has one or more columns, set the values of
            % these defining properties:
            %
            %     FirstCornerLon
            %     DeltaLonNumerator
            %     DeltaLonDenominator
            %
            % to be consistent with a specific LONLIM value, unless
            % the RasterInterpretation is 'postings' and the raster has
            % only one column. The extent in longitude (and intrinsic X) of
            % such a raster is 0, so it requires special handling.
            
            w = lonlim(1);
            e = lonlim(2);
            if self.DeltaLonNumerator > 0
                dlon = e - w;
                self.FirstCornerLon = w;
            else
                dlon = w - e;
                self.FirstCornerLon = e;
            end
            [self.DeltaLonNumerator, self.DeltaLonDenominator] ...
                = simplifyRatio(dlon, diff(self.XLimIntrinsic));
        end
        
        
        function self = setLongitudeProperties(self, ...
                firstCornerLon, deltaLonNumerator, deltaLonDenominator)
            % Set properties controlling longitude vs. intrinsicX
            %
            % Set the following properties as a group:
            %
            %      FirstCornerLon
            %      DeltaLonNumerator
            %      DeltaLonDenominator
            %
            % These properties in combination with RasterSize(2)
            % determine the longitude limits and must be validated
            % together.
            
            % Validate individual inputs
            fname = 'GeoRasterReference.setLongitudeProperties';
            
            validateattributes(firstCornerLon, {'double'}, ...
                {'real','scalar','finite'}, fname, 'firstCornerLon')
            
            validateattributes(deltaLonNumerator, {'double'}, ...
                {'real','scalar','finite','nonzero'}, ...
                fname, 'deltaLonNumerator')
            
            validateattributes(deltaLonDenominator, {'double'}, ...
                {'real','scalar','finite','positive'}, ...
                fname, 'deltaLonDenominator')
            
            % Assign property values
            self.FirstCornerLon = firstCornerLon;
            [self.DeltaLonNumerator, self.DeltaLonDenominator] ...
                = simplifyRatio(deltaLonNumerator, deltaLonDenominator);
        end
        
    end
    
end

%---------------------------- Helper functions ----------------------------

function tf = beyondPole(lat, angleUnits)
% True if lat falls north of 90 degrees N or south of 90
% degrees S. False otherwise, including when lat is
% NaN-valued.
northPole = northPoleLat(angleUnits);
southPole = -northPole;
tf = (lat < southPole | northPole < lat);
end


function lat = northPoleLat(angleUnits)
if strcmp(angleUnits,'degrees')
    lat = 90;
else
    lat = pi/2;
end
end


function lon = wrapToCycle(lon, angleUnits)
% Wrap longitudes to a full cycle, taking into account the
% current unit of angle.
if strcmp(angleUnits,'degrees')
    lon = wrapTo360(lon);
else
    lon = wrapTo2Pi(lon);
end
end


function cycle = fullCycle(angleUnits)
% Return 360 when AngleUnits is degrees.
% Return 2*pi when AngleUnits is radians.
if strcmp(angleUnits,'degrees')
    cycle = 360;
else
    cycle = 2*pi;
end
end
