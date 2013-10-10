%MapRasterReference Reference raster to map coordinates
%
%   A spatialref.MapRasterReference object encapsulates the relationship
%   between a planar map coordinate system and a system of "intrinsic
%   coordinates" anchored to the columns and rows of a 2-D spatially
%   referenced raster grid or image. Typically the raster is sampled
%   regularly in the planar "world X" and "world Y" coordinates of the
%   map system such that the "intrinsic X" and "world X" axes align and
%   likewise with the "intrinsic Y" and "world Y" axes. When this is
%   true the relationship between the two systems is "rectilinear." More
%   generally (and much more rarely), their relationship may be affine,
%   which allows for a possible rotation (and skew). In either case, the
%   sample spacing in from row to row need not equal the sample spacing
%   from column to column. If the raster data set is interpreted as
%   comprising a grid of cells or pixels, the cells or pixels need
%   not be square. In the most general case, they could conceivably be
%   parallelograms, but in practice they are always rectangular.
%
%   MapRasterReference properties:
%      XLimWorld - Limits of raster in world X [xMin xMax]
%      YLimWorld - Limits of raster in world Y [yMin yMax]
%      RasterSize - Number of cells or samples in each spatial dimension
%      RasterInterpretation - Controls handling of raster edges
%      ColumnsStartFrom - Edge where column indexing starts: 'south' or 'north'
%      RowsStartFrom - Edge where row indexing starts: 'west' or 'east'
%
%   MapRasterReference properties (SetAccess = private):
%      DeltaX - Cell width or sample spacing along rows
%      DeltaY - Cell height or sample spacing along columns
%      RasterWidthInWorld - Full extent in along-row direction
%      RasterHeightInWorld - Full extent in along-column direction
%      XLimIntrinsic - Limits of raster in intrinsic X [xMin xMax]
%      YLimIntrinsic - Limits of raster in intrinsic Y [yMin yMax]
%      TransformationType - Transformation type: 'rectilinear' or 'affine'
%      CoordinateSystemType - Type of external system (constant: 'planar')
%
%   MapRasterReference methods:
%      MapRasterReference - Construct spatial.MapRasterReference object
%      sizesMatch - True if object and raster or image are size-compatible
%      intrinsicToWorld - Convert from intrinsic to world coordinates
%      worldToIntrinsic - Convert from world to intrinsic coordinates
%      worldToSub - World coordinates to row and column subscripts
%      contains - True if raster contains points in world coordinate system
%      firstCornerX - World X coordinate of the (1,1) corner of the raster
%      firstCornerY - World Y coordinate of the (1,1) corner of the raster
%      worldFileMatrix - World file parameters for transformation
%
%   See also MAPRASTERREF, spatialref.GeoRasterReference, refmatToMapRasterReference

% Copyright 2010-2011 The MathWorks, Inc.
% $Revision: 1.1.6.4.2.2 $  $Date: 2011/02/05 19:22:34 $

classdef (Sealed = true) MapRasterReference
    
    %------------------- Properties: Public + visible --------------------
    
    properties (Dependent = true)
        %XLimWorld - Limits of raster in world X [xMin xMax]
        %
        %    XLimWorld is a two-element row vector.
        XLimWorld
        
        %YLimWorld - Limits of raster in world Y [yMin yMax]
        %
        %    YLimWorld is a two-element row vector.
        YLimWorld
        
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
        %DeltaX - Cell width or sample spacing along rows
        %
        %   When RasterInterpretation is 'cells', DeltaX equals the
        %   cell width. When RasterInterpretation is 'postings', DeltaX
        %   is the sample spacing along rows (from column to column). In
        %   the case of a rectilinear transformation, DeltaX is signed.
        %   A positive sign indicates that world X increases when
        %   intrinsic X increases, and a negative sign indicates
        %   otherwise. In the case of a general affine transformation
        %   (in which the X-axes may not be aligned), DeltaX also
        %   indicates either cell width or sample spacing along rows,
        %   but it cannot meaningfully be given a sign, so it is
        %   strictly positive.
        DeltaX
        
        %DeltaY - Cell height or sample spacing along columns
        %
        %   When RasterInterpretation is 'cells', DeltaY equals the
        %   cell height. When RasterInterpretation is 'postings', DeltaY
        %   is the sample spacing along columns (from row to row). In
        %   the case of a rectilinear transformation, DeltaY is signed.
        %   A positive sign indicates that world Y increases when
        %   intrinsic Y increases, and a negative sign indicates
        %   otherwise. In the case of a general affine transformation
        %   (in which the Y-axes may not be aligned), DeltaY also
        %   indicates either cell height or sample spacing along columns,
        %   but it cannot meaningfully be given a sign, so it is
        %   strictly positive.
        DeltaY
        
        %RasterWidthInWorld - Full extent in along-row direction
        %
        %   RasterWidthInWorld is the extent of the full raster
        %   or image as measured in the world system in a direction
        %   parallel to its rows. In the case of a rectilinear geometry,
        %   which is most typical, this is the horizontal direction
        %   (east-west).
        RasterWidthInWorld
        
        %RasterHeightInWorld - Full extent in along-column direction
        %
        %   RasterHeightInWorld is the extent of the full raster
        %   or image as measured in the world system in a direction
        %   parallel to its columns. In the case of a rectilinear
        %   geometry, which is most typical, this is the vertical
        %   direction (north-south).
        RasterHeightInWorld
    end
    
    properties (SetAccess = private, Transient = true)
        %XLimIntrinsic - Limits of raster in intrinsic X [xMin xMax]
        %
        %    XLimIntrinsic is a two-element row vector. For an M-by-N
        %    raster with RasterInterpretation equal to 'postings' it
        %    equals [1 N], and for 'cells' it equals [0.5, N + 0.5].
        XLimIntrinsic
        
        %YLimIntrinsic - Limits of raster in intrinsic Y [yMin yMax]
        %
        %    YLimIntrinsic is a two-element row vector. For an M-by-N
        %    raster with RasterInterpretation equal to 'postings' it
        %    equals [1 M], and for 'cells' it equals [0.5, M + 0.5].
        YLimIntrinsic
    end
    
    properties (Dependent = true, SetAccess = private)
        %TransformationType - Transformation type: 'rectilinear' or 'affine'
        %
        %   TransformationType is a string describing the type of geometric
        %   relationship between the intrinsic coordinate system and the
        %   world coordinate system. Its value is 'rectilinear' when world
        %   X depends only on intrinsic X and vice versa, and world Y
        %   depends only on intrinsic Y and vice versa. When the value is
        %   'rectilinear', the image will display without rotation
        %   (although it may be flipped) in the world system. Otherwise the
        %   value is 'affine'.
        TransformationType
    end
    
    properties (Constant = true)
        %CoordinateSystemType - Type of external system (constant: 'planar')
        %
        %   CoordinateSystemType describes the type of coordinate system
        %   to which the image or raster is referenced. It is a constant
        %   string with value 'planar'.
        CoordinateSystemType = 'planar';
    end
    
    %---------------- Properties: Private + hidden ---------------------
    
    properties (Access = private, Hidden = true)
        % The RasterSize and RasterInterpretation properties depend
        % solely on a hidden instance of an IntrinsicRaster2D object. It
        % supports the public sizesMatch(), worldToSub(), and contains()
        % methods, as well as the get methods for XLimWorld, YLimWorld,
        % XCornersWorld, YCornersWorld, and several private methods.
        Intrinsic = spatialref.IntrinsicRaster2D();
        
        % The world limits, column/row direction, delta, raster size, and
        % transformation type properties, along with all the methods except
        % for sizesMatch, depend on a hidden geometric transformation
        % object stored in the Transformation property. It is not
        % initialized, because until the (spatialref.MapRasterReference)
        % constructor runs, we cannot tell if Transformation will hold an
        % instance of a spatialref.RectilinearTransformation object or an
        % instance of a spatialref.AffineTransformation object.
        Transformation
        
        % In a future version, the following will hold the value of a
        % dependent LengthUnits property.
        pLengthUnits = '';
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
        
        function self = MapRasterReference(varargin)
            %MapRasterReference Construct spatialref.MapRasterReference object
            %
            %   R = spatialref.MapRasterReference() constructs a
            %   spatialref.MapRasterReference object with the following
            %   default property settings:
            %
            %                  XLimWorld: [0.5 2.5]
            %                  YLimWorld: [0.5 2.5]
            %                 RasterSize: [2 2]
            %       RasterInterpretation: 'cells'
            %           ColumnsStartFrom: 'south'
            %              RowsStartFrom: 'west'
            %                     DeltaX: 1
            %                     DeltaY: 1
            %         RasterWidthInWorld: 2
            %        RasterHeightInWorld: 2
            %              XLimIntrinsic: [0.5 2.5]
            %              YLimIntrinsic: [0.5 2.5]
            %         TransformationType: 'rectilinear'
            %       CoordinateSystemType: 'planar'
            %
            %   R = spatialref.MapRasterReference(rasterSize, ...
            %       rasterInterpretation, lengthUnits, firstCornerX,...
            %       firstCornerY, jacobianNumerator, jacobianDenominator)
            %   constructs a spatialref.MapRasterReference object from
            %   the following inputs (all 7 must be provided when the
            %   constructor is called):
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
            %     rasterInterpretation -- A string, either 'cells' or
            %        'postings', used to set the RasterInterpretation
            %        property.
            %
            %     lengthUnits -- A place holder to help with forward
            %        compatibility. It must be the empty string, ''.
            %
            %     firstCornerX, firstCornerY -- Scalar values
            %        defining the world X and world Y position of the
            %        outermost corner of the first cell (1,1) of the
            %        raster (if rasterInterpretation is 'cells') or the
            %        first sample point (if rasterInterpretation is
            %        'postings').
            %
            %     jacobianNumerator -- Real, non-singular 2-by-2 matrix
            %     jacobianDenominator -- Real, positive 2-by-2 matrix
            %
            %         The ratio J = jacobianNumerator./jacobianDenominator
            %         determines the Jacobian matrix indicating how
            %         the world coordinates of a point change with
            %         respect to changes in its intrinsic coordinates:
            %
            %         J(1,1) = change in world X wrt intrinsic X
            %         J(1,2) = change in world X wrt intrinsic Y
            %         J(2,1) = change in world Y wrt intrinsic X
            %         J(2,2) = change in world Y wrt intrinsic Y
            %
            %        All 4 elements of jacobianDenominator must be
            %        strictly positive. Because both rectilinear and
            %        affine mappings are both linear, J is invariant
            %        with respect to point location. In the rectilinear
            %        case, J is a diagonal matrix.
            %
            %   See also MAPRASTERREF.
            
            if nargin == 0
                % Default constructor.
                
                % Assign transient limit properties. These two properties
                % are not explicitly dependent, but are architecturally
                % dependent, and don't have default values.
                self.XLimIntrinsic = self.Intrinsic.XLim;
                self.YLimIntrinsic = self.Intrinsic.YLim;
                
                % Construct rectilinear transformation object.
                self.Transformation = spatialref.RectilinearTransformation;
                
                % Set default tie points to be consistent with a
                % single cell, 1-by-1 in world units, centered at
                % xWorld = 1, yWorld = 1, and maintain an identity
                % transformation between the two systems.
                self.Transformation.TiePointIntrinsic = [0.5; 0.5];
                self.Transformation.TiePointWorld     = [0.5; 0.5];
                
                % Set Jacobian matrix to be the 2-by-2 identity matrix.
                self.Transformation.Jacobian = struct( ...
                    'Numerator', [1 0; 0 1], 'Denominator', [1 1; 1 1]);
                
            elseif nargin == 7
                % Construct from complete list of defining properties.
                
                % Deal out values.
                [rasterSize, rasterInterpretation, lengthUnits, ...
                    firstCornerX, firstCornerY, ...
                    jacobianNumerator, jacobianDenominator] = deal(varargin{:});
                
                % Construct (and validate) hidden component object.
                self.Intrinsic = spatialref.IntrinsicRaster2D( ...
                    rasterSize, rasterInterpretation);
                self.Intrinsic.validate()
                
                % Assign transient limit properties.
                self.XLimIntrinsic = self.Intrinsic.XLim;
                self.YLimIntrinsic = self.Intrinsic.YLim;
                
                % Validate corner inputs.
                classname = 'spatialref.MapRasterReference';
                
                % In this version, length units can only be the empty
                % string.
                assert(isequal(lengthUnits,''), ...
                    'map:spatialref:lengthUnitsMustBeEmptyString', ...
                    'Length units must be the empty string, ''''.')
                
                validateattributes(firstCornerX, {'double'}, ...
                    {'real','scalar','finite'}, classname, 'firstCornerX')
                
                validateattributes(firstCornerY, {'double'}, ...
                    {'real','scalar','finite'}, classname, 'firstCornerY')
                
                % Construct transformation object (rectilinear or affine).
                tiePointIntrinsic ...
                    = [self.XLimIntrinsic(1); self.YLimIntrinsic(1)];
                
                tiePointWorld = [firstCornerX; firstCornerY];
                
                self.Transformation = makeGeometricTransformation( ...
                    tiePointIntrinsic, tiePointWorld, ...
                    jacobianNumerator, jacobianDenominator);
            else
                error('map:spatialref:expectedSevenInputs', ...
                    'Either 0 or exactly 7 input arguments are required to construct a %s object.', ...
                    'spatialref.MapRasterReference')
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
        
        
        function [xw, yw] = intrinsicToWorld(self, xi, yi)
            %intrinsicToWorld Convert from intrinsic to world coordinates
            %
            %   [xWorld, yWorld] = R.intrinsicToWorld(...
            %   xIntrinsic, yIntrinsic) maps point locations from the
            %   intrinsic system (xIntrinsic, yIntrinsic) to the world
            %   system (xWorld, yWorld) based on the relationship
            %   defined by the referencing object R. The input may
            %   include values that fall completely outside limits of
            %   the raster (or image) in the intrinsic system. In this
            %   case world X and Y are extrapolated outside the bounds
            %   of the image in the world system.
            
            validateCoordinatePairs(xi, yi, ...
                'spatialref.MapRasterReference.intrinsicToWorld', ...
                'xIntrinsic', 'yIntrinsic')
            
            [xw, yw] = self.Transformation.intrinsicToWorld(xi, yi);
        end
        
        
        function [xi, yi] = worldToIntrinsic(self, xw, yw)
            %worldToIntrinsic Convert from world to intrinsic coordinates
            %
            %   [xIntrinsic, yIntrinsic] = R.worldToIntrinsic(...
            %   xWorld, yWorld) maps point locations from the
            %   world system (xWorld, yWorld) to the intrinsic
            %   system (xIntrinsic, yIntrinsic) based on the relationship
            %   defined by the referencing object R. The input may
            %   include values that fall completely outside limits of
            %   the raster (or image) in the world system. In this
            %   case world X and Y are extrapolated outside the bounds
            %   of the image in the intrinsic system.
            
            validateCoordinatePairs(xw, yw, ...
                'spatialref.MapRasterReference.worldToIntrinsic', ...
                'xWorld', 'yWorld')
            
            [xi, yi] = self.Transformation.worldToIntrinsic(xw, yw);
        end
        
        
        function [row,col] = worldToSub(self, xw, yw)
            %worldToSub World coordinates to row and column subscripts
            %
            %   [I,J] = R.worldToSub(xWorld, yWorld) returns the subscript
            %   arrays I and J. When the referencing object R has
            %   RasterInterpretation 'cells', these are the row and column
            %   subscripts of the raster cells (or image pixels) containing
            %   each element of a set of points given their world
            %   coordinates (xWorld, yWorld).  If R.RasterInterpretation is
            %   'postings', then the subscripts refer to the nearest sample
            %   point (posting). xWorld and yWorld must have the same size.
            %   I and J will have the same size as xWorld and yWorld. For
            %   an M-by-N raster, 1 <= I <= M and 1 <= J <= N, except when
            %   a point xWorld(k), yWorld(k) falls outside the image, as
            %   defined by R.contains(xWorld, yWorld), then
            %   both I(k) and J(k) are NaN.
            
            validateCoordinatePairs(xw, yw, ...
                'spatialref.MapRasterReference.worldToSub', ...
                'xWorld', 'yWorld')
            
            [xi, yi] = self.Transformation.worldToIntrinsic(xw, yw);
            [row, col] = self.Intrinsic.intrinsicToSub(xi, yi);
        end
        
        
        function tf = contains(self, xw, yw)
            %contains True if raster contains points in world coordinate system
            %
            %   TF = R.contains(xWorld, yWorld) returns a logical array TF
            %   having the same size as xWorld, yWorld such that TF(k) is
            %   true if and only if the point (xWorld(k), yWorld(k)) falls
            %   within the bounds of the raster associated with
            %   referencing object R.
            
            validateCoordinatePairs(xw, yw, ...
                'spatialref.MapRasterReference.contains', ...
                'xWorld', 'yWorld')
            
            [xi, yi] = self.Transformation.worldToIntrinsic(xw, yw);
            tf = self.Intrinsic.contains(xi, yi);
        end
        
        
        function xw = firstCornerX(self)
            %firstCornerX - World X coordinate of the (1,1) corner of the raster
            %
            %   R.firstCornerX returns the world X coordinate of the
            %   outermost corner of the first cell (1,1) of the raster
            %   associated with referencing object R (if
            %   R.RasterInterpretation is 'cells') or the first sample
            %   point (if R.RasterInterpretation is 'postings').
            xw = self.Transformation.TiePointWorld(1);
        end
        
        
        function yw = firstCornerY(self)
            %firstCornerY - World Y coordinate of the (1,1) corner of the raster
            %
            %   R.firstCornerY returns the world Y coordinate of the
            %   outermost corner of the first cell (1,1) of the raster
            %   associated with referencing object R (if
            %   R.RasterInterpretation is 'cells') or the first sample
            %   point (if R.RasterInterpretation is 'postings').
            yw = self.Transformation.TiePointWorld(2);
        end
        
        
        function W = worldFileMatrix(self)
            %worldFileMatrix - World file parameters for transformation
            %
            %   W = R.worldFileMatrix returns a 2-by-3 world file matrix.
            %   Each of the 6 elements in W matches one of the lines in a
            %   world file corresponding to the rectilinear or affine
            %   transformation defined by the referencing object R.
            %
            %   Given W with the form:
            %
            %                    W = [A B C;
            %                         D E F],
            %
            %   a point (xi, yi) in intrinsic coordinates maps to a point
            %   (xw, yw) in planar world coordinates like this:
            %
            %         xw = A * (xi - 1) + B * (yi - 1) + C
            %         yw = D * (xi - 1) + E * (yi - 1) + F.
            %
            %   Or, more compactly, [xw yw]' = W * [(xi - 1) (yi - 1) 1]'.
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
            %   The expressions above hold for both affine and rectilinear
            %   transformations, but whenever R.TransformationType is
            %   'rectilinear', B, D, W(2,1) and W(1,2) are identically 0.
            %
            %   See also WORLDFILEREAD, WORLDFILEWRITE.
            
            J = self.Transformation.jacobianMatrix();
            [c, f] = self.intrinsicToWorld(1,1);
            W = [J  [c; f]];
        end
        
    end
    
    %------------------------- Overloaded disp --------------------------
    
    methods
        
        function disp(self)
            % Override the default to display the values of the DeltaX
            % and DeltaY properties as rational numbers.
            
            s = evalc('builtin(''disp'',self)');
            
            if isscalar(self)
                [deltaNumerator, deltaDenominator] ...
                    = self.Transformation.rationalDelta();
                
                if deltaDenominator(1) ~= 1
                    s = replaceValueWithRatio(s, 'DeltaX', ...
                        deltaNumerator(1), deltaDenominator(1));
                end
                
                if deltaDenominator(2) ~= 1
                    s = replaceValueWithRatio(s, 'DeltaY', ...
                        deltaNumerator(2), deltaDenominator(2));
                end
            end
            
            fprintf('%s',s)
        end
        
    end
    
    %-------------------------- Set methods ----------------------------
    
    methods
        
        function self = set.RasterSize(self, rasterSize)
            
            % Current dimensions in intrinsic system.
            currentIntrinsicWidth  = diff(self.XLimIntrinsic);
            currentIntrinsicHeight = diff(self.YLimIntrinsic);
            
            % Update (and validate) RasterSize property.
            self.Intrinsic.RasterSize = rasterSize;
            self.Intrinsic.validate()
            
            % Update cached values (right away).
            self.XLimIntrinsic = self.Intrinsic.XLim;
            self.YLimIntrinsic = self.Intrinsic.YLim;
            
            % Rescale the columns of the Jacobian matrix, as appropriate.
            self = rescaleJacobian(self, ...
                currentIntrinsicWidth, currentIntrinsicHeight);
        end
        
        
        function self = set.RasterInterpretation(self, rasterInterpretation)
            
            % Current dimensions in intrinsic system.
            currentIntrinsicWidth  = diff(self.XLimIntrinsic);
            currentIntrinsicHeight = diff(self.YLimIntrinsic);
            
            % Update (and validate) RasterInterpretation property.
            self.Intrinsic.RasterInterpretation = rasterInterpretation;
            self.Intrinsic.validate()
            
            % Update cached values (right away).
            self.XLimIntrinsic = self.Intrinsic.XLim;
            self.YLimIntrinsic = self.Intrinsic.YLim;
            
            % Rescale the columns of the Jacobian matrix, as appropriate.
            self = rescaleJacobian(self, ...
                currentIntrinsicWidth, currentIntrinsicHeight);
            
            % Reset the intrinsic tie point.
            self.Transformation.TiePointIntrinsic ...
                = [self.XLimIntrinsic(1); self.YLimIntrinsic(1)];
        end
        
        
        function self = set.XLimWorld(self, xLimWorld)           
            validateattributes(xLimWorld, ...
                {'double'}, {'real','row','finite','size', [1 2]}, ...
                'spatialref.MapRasterReference.set.XLimWorld', ...
                'xLimWorld')
            
            assert(xLimWorld(1) < xLimWorld(2), ...
                'map:spatialref:expectedAscendingLimits', ...
                'Elements of %s must be ascending in value.', ...
                'xLimWorld')
            
            currentXLimWorld = self.getXLimWorld();
            
            % Take differences of limits (widths of bounding
            % rectangles); these will be positive numbers.
            difference = diff(xLimWorld);
            currentDifference = diff(currentXLimWorld);
            
            % Scale the first row of the Jacobian matrix to match the
            % change in world X extent.
            J = self.Transformation.Jacobian;
            N = J.Numerator;
            D = J.Denominator;
            N(1,:) = N(1,:) * difference;
            D(1,:) = D(1,:) * currentDifference;
            J.Numerator = N;
            J.Denominator = D;
            self.Transformation.Jacobian = J;
            
            % Reset the X component of the tie point to take care of
            % any translation that is also occurring.
            currentTiePointX = self.Transformation.TiePointWorld(1);
            newTiePointX = xLimWorld(1) ...
                + (currentTiePointX - currentXLimWorld(1)) ...
                * difference / currentDifference;
            
            self.Transformation.TiePointWorld(1) = newTiePointX;
        end
        
        
        function self = set.YLimWorld(self, yLimWorld)
            
            validateattributes(yLimWorld, ...
                {'double'}, {'real','row','finite','size', [1 2]}, ...
                'spatialref.MapRasterReference.set.YLimWorld', ...
                'yLimWorld')
            
            assert(yLimWorld(1) < yLimWorld(2), ...
                'map:spatialref:expectedAscendingLimits', ...
                'Elements of %s must be ascending in value.', ...
                'YLimWorld')
            
            currentYLimWorld = self.getYLimWorld();
            
            % Take differences of limits (widths of bounding
            % rectangle); these will be positive numbers.
            difference = diff(yLimWorld);
            currentDifference = diff(currentYLimWorld);
            
            % Scale the second row of the Jacobian matrix to match the
            % change in world Y extent.
            J = self.Transformation.Jacobian;
            N = J.Numerator;
            D = J.Denominator;
            N(2,:) = N(2,:) * difference;
            D(2,:) = D(2,:) * currentDifference;
            J.Numerator = N;
            J.Denominator = D;
            self.Transformation.Jacobian = J;
            
            % Reset the Y component of the tie point to take care of
            % any translation that is also occurring.
            currentTiePointY = self.Transformation.TiePointWorld(2);
            newTiePointY = yLimWorld(1) ...
                + (currentTiePointY - currentYLimWorld(1)) ...
                * difference / currentDifference;
            
            self.Transformation.TiePointWorld(2) = newTiePointY;
        end
        

        function self = set.ColumnsStartFrom(self, edge)
            edge = validatestring(edge, {'south','north'}, ...
                'spatialref.MapRasterReference.set.ColumnsStartFrom', ...
                'edge');
            
            reverseRasterColumns = xor( ...
                self.columnsStartFromSouth(), strcmp(edge, 'south'));
            if reverseRasterColumns
                % The current (end,1) corner will become the new tie point
                % in the world coordinates.
                [newTiePointX, newTiePointY] ...
                    = self.Transformation.intrinsicToWorld( ...
                          self.XLimIntrinsic(1), self.YLimIntrinsic(2));
                
                self.Transformation.TiePointWorld...
                    = [newTiePointX; newTiePointY];
                
                % Change the sign of the second column of the
                % Jacobian matrix.
                J = self.Transformation.Jacobian;
                J.Numerator(:,2) = -J.Numerator(:,2);
                self.Transformation.Jacobian = J;
            end
        end
        
        
        function self = set.RowsStartFrom(self, edge)
            edge = validatestring(edge, {'east','west'}, ...
                'spatialref.MapRasterReference.set.RowsStartFrom', ...
                'edge');
            
            reverseRasterRows = xor( ...
                self.rowsStartFromWest(), strcmp(edge, 'west'));            
            if reverseRasterRows
                % The current (1,end) corner will become the new tie point
                % in the world coordinates.
                [newTiePointX, newTiePointY] ...
                    = self.Transformation.intrinsicToWorld( ...
                          self.XLimIntrinsic(2), self.YLimIntrinsic(1));
                
                self.Transformation.TiePointWorld...
                    = [newTiePointX; newTiePointY];
                % Change the sign of the first column of the
                % Jacobian matrix.
                J = self.Transformation.Jacobian;
                J.Numerator(:,1) = -J.Numerator(:,1);
                self.Transformation.Jacobian = J;
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
        
        
        function limits = get.XLimWorld(self)
            limits = self.getXLimWorld();
        end
        
        
        function limits = get.YLimWorld(self)
            limits = self.getYLimWorld();
        end
        
        
        function edge = get.ColumnsStartFrom(self)
            if self.columnsStartFromSouth()
                edge = 'south';
            else
                edge = 'north';
            end
        end
        
        
        function edge = get.RowsStartFrom(self)
            if self.rowsStartFromWest()
                edge = 'west';
            else
                edge = 'east';
            end
        end
        
        
        function deltaX = get.DeltaX(self)
            deltaX = self.Transformation.deltaX();
        end
        
        
        function deltaY = get.DeltaY(self)
            deltaY = self.Transformation.deltaY();
        end
        
        
        function width = get.RasterWidthInWorld(self)
            width = abs(self.DeltaX) * diff(self.XLimIntrinsic);
        end
        
        
        function height = get.RasterHeightInWorld(self)
            height = abs(self.DeltaY) * diff(self.YLimIntrinsic);
        end
        
        
        function type = get.TransformationType(self)
            type =  self.Transformation.TransformationType;
        end
        
    end
    
    %---------------------- Private methods ----------------------------
    
    methods (Access = private)
        
        function self = rescaleJacobian(self, ...
                previousIntrinsicWidth, previousIntrinsicHeight)
            % Update the Jacobian matrix in response to a change in the
            % intrinsic dimensions of the raster (which could be due to a
            % change in RasterSize or in RasterInterpretation).
            
            % New dimensions in intrinsic system.
            newIntrinsicWidth  = diff(self.Intrinsic.XLim);
            newIntrinsicHeight = diff(self.Intrinsic.YLim);
            
            % Rescale the columns of the Jacobian matrix, as appropriate.
            J = self.Transformation.Jacobian;
            N = J.Numerator;
            D = J.Denominator;
            
            if previousIntrinsicWidth > 0 && newIntrinsicWidth > 0
                % This is the typical case. Scale the first column of the
                % Jacobian matrix such that the raster continues to fit
                % exactly within the current limits. In all other cases
                % simply leave the Jacobian matrix as-is.
                N(:,1) = N(:,1) * previousIntrinsicWidth;
                D(:,1) = D(:,1) * newIntrinsicWidth;
            end
            
            if previousIntrinsicHeight > 0 && newIntrinsicHeight > 0
                % This is the typical case. Scale the second column of the
                % Jacobian matrix such that the raster continues to fit
                % exactly within the current limits. In all other cases
                % simply leave the Jacobian matrix as-is.
                N(:,2) = N(:,2) * previousIntrinsicHeight;
                D(:,2) = D(:,2) * newIntrinsicHeight;
            end
            
            J.Numerator = N;
            J.Denominator = D;
            self.Transformation.Jacobian = J;
        end
        
        
        function limits = getXLimWorld(self)
            % X-limits of bounding rectangle in world system
            xi = self.XLimIntrinsic([1 1 2 2]);
            yi = self.YLimIntrinsic([1 2 1 2]);
            [xw, ~] = self.Transformation.intrinsicToWorld(xi, yi);
            limits = [min(xw), max(xw)];
        end
        
        
        function limits = getYLimWorld(self)
            % Y-limits of bounding rectangle in world system            
            xi = self.XLimIntrinsic([1 1 2 2]);
            yi = self.YLimIntrinsic([1 2 1 2]);
            [~, yw] = self.intrinsicToWorld(xi, yi);
            limits = [min(yw), max(yw)];
        end
        
        
        function tf = rowsStartFromWest(self)
            % True if and only if rows start from due west,
            % +/- an angle of pi/2.
            
            % Angle between intrinsic X axis and world X axis
            J = self.Transformation.jacobianMatrix();
            alpha = atan2(J(2,1), J(1,1));
            
            tf = (-pi/2 < alpha && alpha <= pi/2);
        end
        
        
        function tf = columnsStartFromSouth(self)
            % True if and only if columns start from due south,
            % +/- an angle of pi/2.
            
            % Angle between intrinsic Y axis and world Y axis
            J = self.Transformation.jacobianMatrix();
            beta = atan2(J(1,2), J(2,2));
            
            tf = -pi/2 < beta && beta <= pi/2;
        end
        
    end
    
end
