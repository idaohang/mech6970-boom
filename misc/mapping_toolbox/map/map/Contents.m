% Mapping Toolbox
% Version 3.3 (R2011a) 09-Mar-2011
%
% GEOSPATIAL DATA IMPORT AND ACCESS
%
% Standard File Formats
%   arcgridread    - Read gridded data set in Arc ASCII Grid Format
%   geotiffinfo    - Information about GeoTIFF file
%   geotiffread    - Read GeoTIFF file
%   geotiffwrite   - Write GeoTIFF file
%   getworldfilename - Derive worldfile name from image file name
%   kmlwrite       - Write geographic data to KML file
%   makeattribspec - Attribute specification from geographic data structure
%   makedbfspec    - DBF specification from geographic data structure
%   sdtsdemread    - Read data from SDTS raster/DEM data set
%   sdtsinfo       - Information about SDTS data set
%   shapeinfo      - Information about shapefile
%   shaperead      - Read vector features and attributes from shapefile
%   shapewrite     - Write geographic data structure to shapefile
%   worldfileread  - Read world file and return referencing object or matrix
%   worldfilewrite - Write world file from referencing object or matrix
%
% Web Map Service
%   wmsfind         - Search local database for Web map servers and layers
%   wmsinfo         - Information about WMS server from capabilities document
%   wmsread         - Retrieve WMS map from server
%   wmsupdate       - Synchronize WMSLayer object with server
%   WebMapServer    - Web map server object
%   WMSCapabilities - Web Map Service capabilities object
%   WMSLayer        - Web Map Service layer object
%   WMSMapRequest   - Web Map Service map request object
%
% Gridded Terrain and Bathymetry Products
%   dted         - Read U.S. Dept. of Defense Digital Terrain Elevation Data (DTED)
%   dteds        - DTED filenames for latitude-longitude quadrangle
%   etopo        - Read gridded global relief data (ETOPO products)
%   globedem     - Read Global Land One-km Base Elevation (GLOBE) data
%   globedems    - GLOBE data filenames for latitude-longitude quadrangle
%   gtopo30      - Read 30-arc-second global digital elevation data (GTOPO30)
%   gtopo30s     - GTOPO30 data filenames for latitude-longitude quadrangle
%   satbath      - Read 2-minute global terrain/bathymetry from Smith and Sandwell
%   tbase        - Read 5-minute global terrain elevations from TerrainBase
%   usgs24kdem   - Read USGS 7.5-minute (30-m or 10-m) Digital Elevation Model
%   usgsdem      - Read USGS 1-degree (3-arc-second) Digital Elevation Model
%   usgsdems     - USGS 1-Degree DEM filenames for latitude-longitude quadrangle
%
% Vector Map Products
%   dcwdata      - Read selected DCW worldwide basemap data
%   dcwgaz       - Search DCW worldwide basemap gazette file
%   dcwread      - Read DCW worldwide basemap file
%   dcwrhead     - Read DCW worldwide basemap file headers
%   fipsname     - Read FIPS name file used with TIGER thinned boundary files
%   gshhs        - Read Global Self-consistent Hierarchical High-resolution Shoreline
%   vmap0data    - Read selected data from Vector Map Level 0
%   vmap0read    - Read Vector Map Level 0 file
%   vmap0rhead   - Read Vector Map Level 0 file headers
%
% Miscellaneous Data Sets
%   avhrrgoode   - Read AVHRR data product stored in Goode projection
%   avhrrlambert - Read AVHRR data product stored in eqaazim projection
%   egm96geoid   - Read 15-minute gridded geoid heights from EGM96
%   readfk5      - Read Fifth Fundamental Catalog of stars
%
%  GUIs for Data Import 
%   demdataui    - UI for selecting digital elevation data
%   vmap0ui      - UI for selecting data from Vector Map Level 0
%
% File Reading Utilities
%   grepfields   - Identify matching records in fixed record length files
%   readfields   - Read fields or records from fixed format file
%   readmtx      - Read matrix stored in file
%   spcread      - Read columns of data from ASCII text file
%
% Ellipsoids, Radii, Areas, and Volumes
%   almanac      - Parameters for Earth, planets, Sun, and Moon
%   earthRadius  - Mean radius of planet Earth
%
% VECTOR MAP DATA AND GEOGRAPHIC DATA STRUCTURES
%
% Geographic Data Structures
%   extractfield - Field values from structure array
%   extractm     - Coordinate data from line or patch display structure
%   updategeostruct - Convert line or patch display structure to geostruct
%
% Data Manipulation
%   bufferm     - Buffer zones for latitude-longitude polygons
%   flatearthpoly - Insert points along dateline to pole
%   interpm     - Densify latitude-longitude sampling in lines or polygons
%   intrplat    - Interpolate latitude at given longitude
%   intrplon    - Interpolate longitude at given latitude
%   ispolycw    - True if polygon vertices are in clockwise order
%   nanclip     - Clip vector data with NaNs at specified pen-down locations
%   poly2ccw    - Convert polygon contour to counterclockwise vertex ordering
%   poly2cw     - Convert polygon contour to clockwise vertex ordering
%   poly2fv     - Convert polygonal region to patch faces and vertices
%   polycut     - Compute branch cuts for holes in polygons
%   polymerge   - Merge line segments with matching endpoints
%   reducem     - Reduce density of points in vector data
%
% Utilities for NaN-Separated Polygons and Lines
%   closePolygonParts  - Close all rings in multipart polygon
%   isShapeMultipart   - True if polygon or line has multiple parts
%   polyjoin    - Convert line or polygon parts from cell arrays to vector form
%   polysplit   - Convert line or polygon parts from vector form to cell arrays
%   removeExtraNanSeparators - Clean up NaN separators in polygons and lines
%
%
% GEOREFERENCED IMAGES AND DATA GRIDS
%
% Spatial Referencing
%   latlon2pix  - Convert latitude-longitude coordinates to pixel coordinates
%   limitm      - Latitude and longitude limits for regular data grid
%   makerefmat  - Construct affine spatial-referencing matrix
%   map2pix     - Convert map coordinates to pixel coordinates
%   mapbbox     - Compute bounding box of georeferenced image or data grid
%   mapoutline  - Compute outline of georeferenced image or data grid
%   meshgrat    - Construct map graticule for surface object display
%   pix2latlon  - Convert pixel coordinates to latitude-longitude coordinates
%   pix2map     - Convert pixel coordinates to map coordinates
%   pixcenters  - Compute pixel centers for georeferenced image or data grid
%   refmatToWorldFileMatrix - Convert referencing matrix to world file matrix
%   refmat2vec  - Convert referencing matrix to referencing vector
%   refvec2mat  - Convert referencing vector to referencing matrix
%   setltln     - Convert data grid rows and columns to latitude-longitude
%   setpostn    - Convert latitude-longitude to data grid rows and columns
%   worldFileMatrixToRefmat - Convert world file matrix to referencing matrix
%
% Spatial Referencing Objects
%   georasterref - Construct spatialref.GeoRasterReference object
%   maprasterref - Construct spatialref.MapRasterReference object
%   refmatToGeoRasterReference - Referencing matrix to GeoRasterReference object
%   refmatToMapRasterReference - Referencing matrix to MapRasterReference object
%   refvecToGeoRasterReference - Referencing vector to GeoRasterReference object
%   spatialref.GeoRasterReference - Reference raster to geographic coordinates
%   spatialref.MapRasterReference - Reference raster to map coordinates
%
% Terrain Analysis
%   gradientm   - Calculate gradient, slope and aspect of data grid
%   los2        - Line of sight visibility between two points in terrain
%   viewshed    - Areas visible from point on terrain elevation grid
%
% Other Analysis/Access
%   areamat     - Surface area covered by non-zero values in binary data grid
%   filterm     - Filter latitudes/longitudes based on underlying data grid
%   findm       - Latitudes and longitudes of non-zero data grid elements
%   ltln2val    - Extract data grid values for specified locations
%   mapprofile  - Interpolate between waypoints on regular data grid
%
% Construction and Modification
%   changem     - Substitute values in data array
%   encodem     - Fill in regular data grid from seed values and locations
%   geoloc2grid - Convert geolocated data array to regular data grid
%   imbedm      - Encode data points into regular data grid
%   neworig     - Orient regular data grid to oblique aspect
%   resizem     - Resize regular data grid
%   sizem       - Row and column dimensions needed for regular data grid
%   vec2mtx     - Convert latitude-longitude vectors to regular data grid
%
% Initialization
%   nanm        - Construct regular data grid of NaNs
%   onem        - Construct regular data grid of 1s
%   spzerom     - Construct sparse regular data grid of 0s
%   zerom       - Construct regular data grid of 0s
%
%
% MAP PROJECTIONS AND COORDINATES
%
% Available Map Projections (in addition to PROJ.4 library)
%   maps        - List available map projections and verify names
%   maplist     - Map projections available in Mapping Toolbox
%   projlist    - Map projections supported by PROJFWD and PROJINV
%
% Map Projection Transformations
%   mfwdtran    - Project geographic features to map coordinates
%   minvtran    - Unproject features from map to geographic coordinates
%   projfwd     - Forward map projection using PROJ.4 library
%   projinv     - Inverse map projection using PROJ.4 library
%
% Angles, Scales and Distortions
%   vfwdtran    - Direction angle in map plane from azimuth on ellipsoid
%   vinvtran    - Azimuth on ellipsoid from direction angle in map plane
%   distortcalc - Distortion parameters for map projections
%
% Visualizing Map Distortions
%   mdistort    - Display contours of constant map distortion
%   tissot      - Project Tissot indicatrices on map axes
%
% Cylindrical Projections
%   balthsrt    - Balthasart Cylindrical Projection
%   behrmann    - Behrmann Cylindrical Projection
%   bsam        - Bolshoi Sovietskii Atlas Mira Cylindrical Projection
%   braun       - Braun Perspective Cylindrical Projection
%   cassini     - Cassini Transverse Cylindrical Projection
%   cassinistd  - Cassini Transverse Cylindrical Projection -- Standard
%   ccylin      - Central Cylindrical Projection
%   eqacylin    - Equal Area Cylindrical Projection
%   eqdcylin    - Equidistant Cylindrical Projection
%   giso        - Gall Isographic Cylindrical Projection
%   gortho      - Gall Orthographic Cylindrical Projection
%   gstereo     - Gall Stereographic Cylindrical Projection
%   lambcyln    - Lambert Equal Area Cylindrical Projection
%   mercator    - Mercator Cylindrical Projection
%   miller      - Miller Cylindrical Projection
%   pcarree     - Plate Carree Cylindrical Projection
%   tranmerc    - Transverse Mercator Projection
%   trystan     - Trystan Edwards Cylindrical Projection
%   wetch       - Wetch Cylindrical Projection
%
% Pseudocylindrical Projections
%   apianus     - Apianus II Pseudocylindrical Projection
%   collig      - Collignon Pseudocylindrical Projection
%   craster     - Craster Parabolic Pseudocylindrical Projection
%   eckert1     - Eckert I Pseudocylindrical Projection
%   eckert2     - Eckert II Pseudocylindrical Projection
%   eckert3     - Eckert III Pseudocylindrical Projection
%   eckert4     - Eckert IV Pseudocylindrical Projection
%   eckert5     - Eckert V Pseudocylindrical Projection
%   eckert6     - Eckert VI Pseudocylindrical Projection
%   flatplrp    - McBryde-Thomas Flat-Polar Parabolic Projection
%   flatplrq    - McBryde-Thomas Flat-Polar Quartic Projection
%   flatplrs    - McBryde-Thomas Flat-Polar Sinusoidal Projection
%   fournier    - Fournier Pseudocylindrical Projection
%   goode       - Goode Homolosine Pseudocylindrical Projection
%   hatano      - Hatano Asymmetrical Equal Area Pseudocylindrical Projection
%   kavrsky5    - Kavraisky V Pseudocylindrical Projection
%   kavrsky6    - Kavraisky VI Pseudocylindrical Projection
%   loximuth    - Loximuthal Pseudocylindrical Projection
%   modsine     - Tissot Modified Sinusoidal Pseudocylindrical Projection
%   mollweid    - Mollweide Pseudocylindrical Projection
%   putnins5    - Putnins P5 Pseudocylindrical Projection
%   quartic     - Quartic Authalic Pseudocylindrical Projection
%   robinson    - Robinson Pseudocylindrical Projection
%   sinusoid    - Sinusoidal Pseudocylindrical Projection
%   wagner4     - Wagner IV Pseudocylindrical Projection
%   winkel      - Winkel I Pseudocylindrical Projection
%
% Conic Projections
%   eqaconic    - Albers Equal Area Conic Projection
%   eqaconicstd - Albers Equal Area Conic Projection -- Standard
%   eqdconic    - Equidistant Conic Projection
%   eqdconicstd - Equidistant Conic Projection -- Standard
%   lambert     - Lambert Conformal Conic Projection
%   lambertstd  - Lambert Conformal Conic Projection -- Standard
%   murdoch1    - Murdoch I Conic Projection
%   murdoch3    - Murdoch III Minimum Error Conic Projection
%
% Polyconic and Pseudoconic Projections
%   bonne       - Bonne Pseudoconic Projection
%   polycon     - Polyconic Projection
%   polyconstd  - Polyconic Projection -- Standard
%   vgrint1     - Van Der Grinten I Polyconic Projection
%   werner      - Werner Pseudoconic Projection
%
% Azimuthal, Pseudoazimuthal and Modified Azimuthal Projections
%   aitoff      - Aitoff Modified Azimuthal Projection
%   breusing    - Breusing Harmonic Mean Azimuthal Projection
%   bries       - Briesemeister's Modified Azimuthal Projection
%   eqaazim     - Lambert Equal Area Azimuthal Projection
%   eqdazim     - Equidistant Azimuthal Projection
%   gnomonic    - Gnomonic Azimuthal Projection
%   hammer      - Hammer Modified Azimuthal Projection
%   ortho       - Orthographic Azimuthal Projection
%   stereo      - Stereographic Azimuthal Projection
%   vperspec    - Vertical Perspective Azimuthal Projection
%   wiechel     - Wiechel Equal Area Pseudoazimuthal Projection
%
% UTM and UPS Systems
%   ups         - Universal Polar Stereographic system
%   utm         - Universal Transverse Mercator system
%   utmgeoid    - Select ellipsoid for given UTM zone
%   utmzone     - Select UTM zone given latitude and longitude
%
% Three-Dimensional Globe Display
%   globe       - Earth as sphere in 3-D graphics
%
% Rotating Coordinates on the Sphere
%   newpole     - Origin vector to place specific point at pole
%   org2pol     - Location of north pole on rotated map
%   putpole     - Origin vector to place north pole at specific point
%
% Map Trimming
%   maptriml    - Trim lines to latitude-longitude quadrangle
%   maptrimp    - Trim polygons to latitude-longitude quadrangle
%   maptrims    - Trim regular data grid to latitude-longitude quadrangle
%
% Trimming and Clipping
%   clipdata    - Clip data to +/- pi in longitude, +/- pi/2 in latitude
%   trimcart    - Trim graphic objects to map frame
%   trimdata    - Trim map data exceeding projection limits
%   undoclip    - Remove object clips introduced by CLIPDATA
%   undotrim    - Remove object trims introduced by TRIMDATA
%
%
% MAP DISPLAY AND INTERACTION
%
% Map Creation and High-Level Display
%   axesm       - Define map axes and set map properties
%   displaym    - Display geographic data from display structure
%   geoshow     - Display map latitude and longitude data
%   grid2image  - Display regular data grid as image
%   mapshow     - Display map data without projection
%   mapview     - Interactive map viewer
%   usamap      - Construct map axes for United States of America
%   worldmap    - Construct map axes for given region of world
%
% Vector Symbolization
%   makesymbolspec  - Construct vector symbolization specification 
%
% Displaying Lines and Contours
%   contourm    - Project 2-D contour plot of map data
%   contour3m   - Project 3-D contour plot of map data
%   contourfm   - Project filled 2-D contour plot of map data 
%   linem       - Project line object on map axes
%   plotm       - Project 2-D lines and points on map axes
%   plot3m      - Project 3-D lines and points on map axes
%
% Displaying Patch Data
%   fillm       - Project filled 2-D patch objects on map axes
%   fill3m      - Project filled 3-D patch objects on map axes
%   patchesm    - Project patches on map axes as individual objects
%   patchm      - Project patch objects on map axes
%
% Displaying Data Grids
%   meshm       - Project regular data grid on map axes
%   pcolorm     - Project regular data grid on map axes in z = 0 plane
%   surfacem    - Project and add geolocated data grid to current map axes
%   surfm       - Project geolocated data grid on map axes
%
% Displaying Light Objects and Lighted Surfaces
%   lightm      - Project light objects on map axes
%   meshlsrm    - 3-D lighted shaded relief of regular data grid
%   surflm      - 3-D shaded surface with lighting on map axes
%   surflsrm    - 3-D lighted shaded relief of geolocated data grid
%   shaderel    - Construct cdata and colormap for shaded relief
%
% Displaying Thematic Maps
%   cometm      - Project 2-D comet plot on map axes
%   comet3m     - Project 3-D comet plot on map axes
%   quiverm     - Project 2-D quiver plot on map axes
%   quiver3m    - Project 3-D quiver plot on map axes
%   scatterm    - Project point markers with variable color and area
%   stem3m      - Project stem plot on map axes
%   symbolm     - Project point markers with variable size
%
% Annotating Map Displays
%   clabelm     - Add contour labels to map contour display
%   clegendm    - Add legend labels to map contour display
%   framem      - Toggle and control display of map frame
%   gridm       - Toggle and control display of graticule lines
%   lcolorbar   - Append colorbar with text labels
%   mlabel      - Toggle and control display of meridian labels
%   mlabelzero22pi - Convert meridian labels to 0 to 360-degree range
%   northarrow  - Add graphic element pointing to geographic North Pole
%   plabel      - Toggle and control display of parallel labels
%   rotatetext  - Rotate text to projected graticule
%   scaleruler  - Add or modify graphic scale on map axes
%   textm       - Project text annotation on map axes
%
% Colormaps for Map Displays
%   contourcmap - Contour colormap and colorbar for current axes
%   demcmap     - Colormaps appropriate to terrain elevation data
%   polcmap     - Colormaps appropriate to political regions
%
% Interactive Map Positions
%   gcpmap      - Get current mouse point from map axes
%   gtextm      - Place text on map using mouse
%   inputm      - Return latitudes and longitudes of mouse click locations
%
% Interactive Track and Circle Definition
%   scircleg    - Small circle defined via mouse input
%   sectorg     - Sector of small circle defined via mouse input
%   trackg      - Great circle or rhumb line defined via mouse input
%
% Graphical User Interfaces
%   axesmui     - Interactively define map axes properties
%   clrmenu     - Add colormap menu to figure window
%   colorm      - Create index map colormaps
%   getseeds    - Interactively assign seeds for data grid encoding
%   lightmui    - Control position of lights on globe or 3-D map
%   maptrim     - Customize map data sets
%   maptool     - Add menu activated tools to map figure
%   mlayers     - Control plotting of display structure elements
%   mobjects    - Manipulate object sets displayed on map axes
%   originui    - Interactively modify map origin
%   panzoom     - Pan and zoom on map axes
%   parallelui  - Interactively modify map parallels
%   qrydata     - Create queries associated with map axes
%   rootlayr    - Construct cell array of workspace variables for MLAYERS tool
%   seedm       - Interactively fill regular data grids with seed values
%   surfdist    - Interactive distance, azimuth and reckoning calculations
%   uimaptbx    - Process button down callbacks for mapped objects
%   utmzoneui   - Choose or identify UTM zone by clicking on map
%
% Map Object and Projection Properties
%   cart2grn    - Transform projected coordinates to Greenwich system
%   defaultm    - Initialize or reset map projection structure
%   gcm         - Current map projection structure
%   geotiff2mstruct - Convert GeoTIFF information to map projection structure
%   getm        - Map object properties
%   handlem     - Handles of displayed map objects
%   ismap       - True for axes with map projection
%   ismapped    - True if object is projected on map axes
%   makemapped  - Convert ordinary graphics object to mapped object
%   namem       - Determine names for valid map graphics objects
%   project     - Project displayed map graphics object
%   restack     - Restack objects within map axes
%   rotatem     - Transform map data to new origin and orientation
%   setm        - Set properties of map axes and graphics objects
%   tagm        - Set tag property of map graphics objects
%   zdatam      - Adjust z-plane of displayed map objects
%
% Controlling Map Appearance
%   axesscale   - Resize axes for equivalent scale
%   camposm     - Set camera position using geographic coordinates
%   camtargm    - Set camera target using geographic coordinates
%   camupm      - Set camera up vector using geographic coordinates
%   daspectm    - Control vertical exaggeration in map display
%   paperscale  - Set figure properties for printing at specified map scale
%   previewmap  - Preview map at printed size
%   tightmap    - Remove white space around map
%
% Clearing Map Displays/Managing Visibility
%   clma        - Clear current map axes
%   clmo        - Clear specified graphic objects from map axes
%   hidem       - Hide specified graphic objects on map axes
%   showaxes    - Toggle display of map coordinate axes
%   showm       - Show specified graphic objects on map axes
%
%
% GEOGRAPHIC CALCULATIONS
%
% Geometry of Sphere and Ellipsoid
%   antipode    - Point on opposite side of globe
%   areaint     - Surface area of polygon on sphere or ellipsoid
%   areaquad    - Surface area of latitude-longitude quadrangle
%   azimuth     - Azimuth between points on sphere or ellipsoid
%   departure   - Departure of longitudes at specific latitudes
%   distance    - Distance between points on sphere or ellipsoid
%   ellipse1    - Geographic ellipse from center, semimajor axes, eccentricity and azimuth
%   gc2sc       - Center and radius of great circle
%   meridianarc - Ellipsoidal distance along meridian
%   meridianfwd - Reckon position along meridian
%   reckon      - Point at specified azimuth, range on sphere or ellipsoid
%   scircle1    - Small circles from center, range and azimuth
%   scircle2    - Small circles from center and perimeter
%   track1      - Geographic tracks from starting point, azimuth and range
%   track2      - Geographic tracks from starting and ending points
%
% Three-Dimensional Coordinates
%   ecef2geodetic - Convert geocentric (ECEF) to geodetic coordinates
%   ecef2lv       - Convert geocentric (ECEF) to local vertical coordinates
%   elevation     - Local vertical elevation angle, range, and azimuth
%   geodetic2ecef - Convert geodetic to geocentric (ECEF) coordinates
%   lv2ecef       - Convert local vertical to geocentric (ECEF) coordinates
%
% Ellipsoids and Latitudes
%   axes2ecc    - Eccentricity of ellipse with given axes lengths
%   convertlat  - Convert between geodetic and auxiliary latitudes
%   ecc2flat    - Flattening of ellipse with given eccentricity
%   ecc2n       - n-value of ellipse with given eccentricity
%   flat2ecc    - Eccentricity of ellipse with given flattening
%   geocentric2geodeticLat - Convert geocentric to geodetic latitude
%   geodetic2geocentricLat - Convert geodetic to geocentric latitude
%   majaxis     - Semimajor axis of ellipse with given semiminor axis and eccentricity
%   minaxis     - Semiminor axis of ellipse with given semimajor axis and eccentricity
%   n2ecc       - Eccentricity of ellipse with given n-value
%   rcurve      - Radii of curvature of ellipsoid
%   rsphere     - Radii of auxiliary spheres
%
% Overlaying Geometric Objects
%   circcirc    - Intersections of circles in Cartesian plane
%   gcxgc       - Intersection points for pairs of great circles
%   gcxsc       - Intersection points for great and small circle pairs
%   ingeoquad        - True for points inside or on lat-lon quadrangle
%   intersectgeoquad - Intersection of two latitude-longitude quadrangles
%   linecirc    - Intersections of circles and lines in Cartesian plane
%   outlinegeoquad   - Polygon outlining geographic quadrangle
%   polybool    - Set operations on polygonal regions
%   polyxpoly   - Intersection points for lines or polygon edges
%   rhxrh       - Intersection points for pairs of rhumb lines
%   scxsc       - Intersection points for pairs of small circles
%
% Geographic Statistics
%   combntns    - All possible combinations of a set of values
%   eqa2grn     - Convert from equal area to Greenwich coordinates
%   grn2eqa     - Convert from Greenwich to equal area coordinates
%   hista       - Histogram for geographic points with equal-area bins
%   histr       - Histogram for geographic points with equirectangular bins
%   meanm       - Mean location of geographic points
%   stdist      - Standard distance for geographic points
%   stdm        - Standard deviation for geographic points
%
% Navigation
%   crossfix    - Cross fix positions from bearings and ranges
%   dreckon     - Dead reckoning positions for track
%   driftcorr   - Heading to correct for wind or current drift
%   driftvel    - Wind or current velocity from heading, course, and speeds
%   gcwaypts    - Equally spaced waypoints along great circle track
%   legs        - Courses and distances between navigational waypoints
%   navfix      - Mercator-based navigational fix
%   timezone    - Time zone based on longitude
%   track       - Track segments to connect navigational waypoints
%
%
% UTILITIES
%
% Angle Conversions
%   degtorad     - Convert angles from degrees to radians
%   degrees2dm   - Convert degrees to degrees-minutes
%   degrees2dms  - Convert degrees to degrees-minutes-seconds
%   dm2degrees   - Convert degrees-minutes to degrees
%   dms2degrees  - Convert degrees-minutes-seconds to degrees
%   fromDegrees  - Convert angles from degrees
%   fromRadians  - Convert angles from radians
%   radtodeg     - Convert angles from radians to degrees
%   str2angle    - Convert strings to angles in degrees
%   toDegrees    - Convert angles to degrees
%   toRadians    - Convert angles to radians
%
% Distance Conversions
%   deg2km       - Convert distance from degrees to kilometers
%   deg2nm       - Convert distance from degrees to nautical miles
%   deg2sm       - Convert distance from degrees to statute miles
%   km2deg       - Convert distance from kilometers to degrees
%   km2nm        - Convert distance from kilometers to nautical miles
%   km2rad       - Convert distance from kilometers to radians
%   km2sm        - Convert distance from kilometers to statute miles
%   nm2deg       - Convert distance from nautical miles to degrees
%   nm2km        - Convert distance from nautical miles to kilometers
%   nm2rad       - Convert distance from nautical miles to radians
%   nm2sm        - Convert distance from nautical to statute miles
%   rad2km       - Convert distance from radians to kilometers
%   rad2nm       - Convert distance from radians to nautical miles
%   rad2sm       - Convert distance from radians to statute miles
%   sm2deg       - Convert distance from statute miles to degrees
%   sm2km        - Convert distance from statute miles to kilometers
%   sm2nm        - Convert distance from statute to nautical miles
%   sm2rad       - Convert distance from statute miles to radians
%
% Conversion Factors for Angles and Distances
%   unitsratio   - Unit conversion factors
%
% Data Precision
%   epsm        - Accuracy in angle units for certain map computations
%   roundn      - Round to multiple of 10^n
%
% Image Conversion
%   ind2rgb8    - Convert indexed image to uint8 RGB image
%
% Longitude or Azimuth Wrapping
%  unwrapMultipart - Unwrap vector of angles with NaN-delimited parts
%  wrapTo180    - Wrap angle in degrees to [-180 180]
%  wrapTo360    - Wrap angle in degrees to [0 360]
%  wrapToPi     - Wrap angle in radians to [-pi pi]  
%  wrapTo2Pi    - Wrap angle in radians to [0 2*pi]
%
% String Formatters
%   angl2str     - Format angle strings
%   dist2str     - Format distance strings
%
% See also MAPDEMOS.

% Copyright 1996-2011 The MathWorks, Inc.
% Generated from Contents.m_template revision 1.1.4.37.2.1  $Date: 2010/12/03 21:43:39 $

% Deprecated and obsolete functions: MAP
%   angledim    - Convert angle units
%   deg2rad     - Convert angles from degrees to radians
%   distdim     - Convert distance units
%   eastof      - Wrap longitudes to values east of specified meridian
%   npi2pi      - Wrap latitudes to [-180 180] degree interval
%   rad2deg     - Convert angles from radians to degrees
%   smoothlong  - Remove discontinuities in longitude data
%   unitstr     - Check unit strings or abbreviations
%   westof      - Wrap longitudes to values west of specified meridian
%   zero22pi    - Wrap longitudes to [0 360) degree interval
%
% Deprecated and obsolete functions: MAPDISP
%   colorui     - Interactively define RGB color
%
% Deprecated and obsolete functions: MAPFORMATS
%   etopo5      - Read 5-minute gridded terrain/bathymetry from global ETOPO5 data set
%   tgrline     - Read TIGER/Line data
%
% Undocumented functions: MAP
%   elpcalc     - Volume and surface area of an oblate spheroid
%   eqacalc     - Transform data to/from an equal area space
%   merccalc    - Transform data to/from a Mercator space
%   sphcalc     - Compute volume and surface area for a sphere
%   mapgate     - Gateway routine to call private functions
%   geoidtst         - Test for a valid geoid vector
%   ignoreComplex    - Convert complex input to real and issue warning
%   num2ordinal      - Convert positive integer to ordinal string
%
% Undocumented functions: MAPDISP (used by MAPTOOL)
%   scirclui    - Interactive tool for adding small circles to a map
%   trackui     - Interactive tool for adding great circles and rhumb lines to a map
%   maphlp1     - Help Utility for Selected GUIs
%   maphlp2     - Help Utility for Selected GUIs
%   maphlp3     - Help Utility for Selected GUIs
%   maphlp4     - Help Utility for Selected GUIs
%   clrpopup    - Processes callback from color popup menus
%   varpick     - Modal pick list to select a variable from the workspace
%
% Undocumented functions: MAPDISP (used outside MAPDISP)
%   degchar     - Return the LaTeX degree symbol character
%   leadblnk    - Delete leading characters common to all rows of a string matrix
%   shiftspc    - Left or right justify a string matrix
%
% Undocumented functions: MAPUTILS
%   checkangleunits - Check and standardize angle units string
%   checkellipsoid  - Check validity of reference ellipsoid vector
%   checkgeoquad    - Validate limits of geographic quadrangle
%   checklatlon     - Validate pair of latitude-longitude arrays
%   checkrefmat     - Check validity of referencing matrix
%   checkrefvec     - Check validity of referencing vector
%
% Deleted functions: MAP
%   deg2dm      - Convert angles from degrees to deg:min encoding
%   deg2dms     - Convert angles from degrees to deg:min:sec encoding
%   dms2deg     - Convert angles from deg:min:sec encoding to degrees
%   dms2dm      - Convert angles from deg:min:sec to deg:min encoding
%   dms2mat     - Expand deg:min:sec encoded vector to [deg min sec] matrix
%   dms2rad     - Convert angles from deg:min:sec encoding to radians
%   hms2hm      - Convert time from hrs:min:sec to hrs:min encoding
%   hms2hr      - Convert time from hrs:min:sec encoding to hours
%   hms2mat     - Expand hrs:min:sec encoded vector to [hrs min sec] matrix
%   hms2sec     - Convert time from hrs:min:sec encoding to seconds
%   hr2hm       - Convert time from hours to hrs:min encoding
%   hr2hms      - Convert time from hours to hrs:min:sec encoding
%   hr2sec      - Convert time from hours to seconds
%   mat2dms     - Collapse [deg min sec] matrix to deg:min:sec encoding
%   mat2hms     - Collapse [hrs min sec] matrix to hrs:min:sec encoding
%   rad2dm      - Convert angles from radians to deg:min encoding
%   rad2dms     - Convert angles from radians to deg:min:sec encoding
%   sec2hm      - Convert time from seconds to hrs:min encoding
%   sec2hms     - Convert time from seconds to hrs:min:sec encoding
%   sec2hr      - Convert time from seconds to hours
%   time2str    - Format time strings
%   timedim     - Convert time units or encodings
%
% Deleted functions: MAPDISP
%   contorm     - Project a contour plot of data onto the current map axes
%   contor3m    - Project a 3D contour plot of data onto the current map axes
%   cmapui      - Create custom colormap
%
% Deleted functions: MAPFORMATS
%   tigermif    - Read a TIGER MIF thinned boundary file
%   tigerp      - Read TIGER p and pa thinned boundary files
