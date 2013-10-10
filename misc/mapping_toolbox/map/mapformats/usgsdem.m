function [map,maplegend,Astruc] = usgsdem(varargin)
%USGSDEM  Read USGS 1-degree (3-arc-second) Digital Elevation Model
%
%  [Z, REFVEC] = USGSDEM(FILENAME, SAMPLEFACTOR) reads the specified file
%  and returns the elevation data in the regular data grid, Z, along with
%  referencing vector REFVEC.  REFVEC is a 1-by-3 vector containing
%  elements [cells/degree north-latitude west-longitude] with latitude and
%  longitude limits in degrees.  The data can be read at full resolution
%  (SAMPLEFACTOR = 1), or can be downsampled by SAMPLEFACTOR. A
%  SAMPLEFACTOR of 3 returns every third point, for example, giving 1/3 of
%  the full resolution.  The grid for the digital elevation maps is based
%  on World Geodetic System 1984 (WGS84). Older DEMs were based on WGS72.
%
%  [Z, REFVEC] = USGSDEM(FILENAME, SAMPLEFACTOR, LATLIM, LONLIM) reads a
%  subset of the elevation data from FILENAME. The limits of the desired
%  data are specified as two element vectors of latitude, LATLIM, and
%  longitude, LONLIM, in degrees. The elements of LATLIM and LONLIM must be
%  in ascending order.  The data will extend somewhat outside the requested
%  area. If omitted, the entire area covered by the DEM file is returned.
%
%  [Z, REFVEC, HEADER] = USGSDEM(...) returns file header in a structure,
%  HEADER.
%
%  For details on locating USGSDEM data for download over the Internet, see
%  the following documentation at the MathWorks web site:
%
%  <a href="matlab:
%  web('http://www.mathworks.com/support/tech-notes/2100/2101.html#demopt')
%  ">http://www.mathworks.com/support/tech-notes/2100/2101.html</a>
%
%  See also DEMDATAUI, USGSDEMS, USGS24KDEM.

% Copyright 1996-2009 The MathWorks, Inc.
% $Revision: 1.1.6.8 $ $Date: 2009/03/30 23:39:28 $
% Written by:  A. Kim

%  Ascii data file
%  Data arranged in S-N rows by W-E columns
%  Elevation in meters

[map,maplegend,Astruc] = usgsdemf(varargin{:});

%--------------------------------------------------------------------------

function [map,maplegend,Astruc] = usgsdemf(fname,scalefactor,latlim,lonlim)

if nargin==2
    subset = 0;
elseif nargin==4
    subset = 1;
else
    error(['map:' mfilename ':invalidArgCount'], ...
        'Incorrect number of arguments')
end

sf = scalefactor;
arcsec3 = 3/60^2;
celldim = sf*arcsec3;
halfcell = celldim/2;

fid = fopen(fname,'r');
if fid==-1
    [fname, path] = uigetfile('', 'select the USGS-DEM file');
    if fname == 0
        return
    end
    fname = [path fname];
    fid = fopen(fname,'r');
end

%  --- Read Record Type A (Header Info) ---

% Define the structure of the header records and data records
Afield = Adescription;

% Read the header record. Some fields may be empty, so assume
% a fixed format as documented in the DEM Data User's Guide 5 (1993),
% and read the header field-by-field.
Astruc = readfields(fname,Afield,1,'native',fid);

% Extract needed data from the header structure
if Astruc.PlanimetricReferenceSystemCode==1
    fclose(fid);
    error(['map:' mfilename ':invalidFileFormat'], ...
        'Ground planimetric coordinates not in arc-seconds. Try USGS24KDEM')
elseif Astruc.PlanimetricReferenceSystemCode~=0
    fclose(fid);
    error(['map:' mfilename ':invalidPlanimetricCoordinates'], ...
        'Ground planimetric coordinates not in arc-seconds.')
end

corners = ( reshape(Astruc.BoundingBox,[2 4])' )/60^2;
ncols = Astruc.NrowsCols(2);
if ~subset
    % Check to see if ncols fit scalefactor
    if mod((ncols-1),sf)~=0
        fclose(fid);
        error(['map:' mfilename ':sampleFactorNotDivisibleIntoCols'], ...
            'SAMPLEFACTOR does not divide evenly into %s columns', ...
            num2str(ncols))
    end
end

dy = arcsec3;
switch ncols
    case 1201, dx = arcsec3;
    case 601,  dx = 2*arcsec3;
    case 401,  dx = 3*arcsec3;
    otherwise
        fclose(fid);
        error(['map:' mfilename ':invalidNcols'], 'Invalid ncols')
end

%  Define border of map
maplatlim(1) = corners(1,2) - halfcell;
maplatlim(2) = corners(2,2) + halfcell;
maplonlim(1) = corners(1,1) - halfcell;
maplonlim(2) = corners(4,1) + halfcell;

if subset

    % Check to see if latlim and lonlim within map limits
    if latlim(1) > latlim(2)
        fclose(fid);
        error(['map:' mfilename ':latlimNotAscending'], ...
            'First element of latlim must be less than second')
    end
    if lonlim(1) > lonlim(2)
        fclose(fid);
        error(['map:' mfilename ':lonlimNotAscending'], ...
            'First element of lonlim must be less than second')
    end

    if  latlim(1) > maplatlim(2) || latlim(2) < maplatlim(1) || ...
        lonlim(1) > maplonlim(2) || lonlim(2) < maplonlim(1)
        warning(['map:' mfilename ':limitsExcludeDataset'], ...
            ['Requested latitude or longitude limits are off the map\n' ...
            ' latlim for this dataset is [%.4f %.4f] \n' ...
            ' lonlim for this dataset is [%.4f %.4f]'], ...
            maplatlim(1), maplatlim(2), maplonlim(1), maplonlim(2))
        map=[];
        maplegend = [];
        fclose(fid);
        return
    end

    clampLimits = false;
    if latlim(1) < maplatlim(1)
        latlim(1) = maplatlim(1);
        clampLimits = true;
    end
    if latlim(2) > maplatlim(2)
        latlim(2) = maplatlim(2);
        clampLimits = true;
    end
    if lonlim(1) < maplonlim(1)
        lonlim(1) = maplonlim(1);
        clampLimits = true;
    end
    if lonlim(2) > maplonlim(2)
        lonlim(2) = maplonlim(2);
        clampLimits = true;
    end
    if clampLimits
        warning(['map:' mfilename ':clampingLimits'], ...
            ['Requested latitude or longitude limits exceed map limits\n' ...
            ' latlim for this dataset is [%.4f %.4f] \n' ...
            ' lonlim for this dataset is [%.4f %.4f]'], ...
            maplatlim(1), maplatlim(2), maplonlim(1), maplonlim(2))
    end

    % Convert lat and lon limits to row and col limits
    halfdy = dy/2;
    halfdx = dx/2;
    ltlwr = corners(1,2)-halfdy:dy:corners(2,2)-halfdy;
    ltupr = corners(1,2)+halfdy:dy:corners(2,2)+halfdy;
    lnlwr = corners(1,1)-halfdx:dx:corners(4,1)-halfdx;
    lnupr = corners(1,1)+halfdx:dx:corners(4,1)+halfdx;
    if latlim(1)>=maplatlim(1) && latlim(1)<=ltlwr(1)
        rowlim(1) = 1;
    else
        rowlim(1) = find(ltlwr<=latlim(1) & ltupr>=latlim(1), 1 );
    end
    if latlim(2)<=maplatlim(2) && latlim(2)>=ltupr(length(ltupr))
        rowlim(2) = 1201;
    else
        rowlim(2) = find(ltlwr<=latlim(2) & ltupr>=latlim(2), 1, 'last' );
    end
    if lonlim(1)==maplonlim(1)
        collim(1) = 1;
    else
        collim(1) = find(lnlwr<=lonlim(1) & lnupr>=lonlim(1), 1 );
    end
    if lonlim(2)==maplonlim(2)
        collim(2) = ncols;
    else
        collim(2) = find(lnlwr<=lonlim(2) & lnupr>=lonlim(2), 1, 'last' );
    end

end

%  --- Read Record Type B (Elevation Data) ---

% Start profile position indicators
startprofiles = 1029:8192:1029+(ncols-1)*8192;
sfmin = 1200/(ncols-1);
if mod(sf,sfmin)~=0
    fclose(fid);
    error(['map:' mfilename ':invalidSamplefactor'], ...
        'Samplefactor must be multiple of %s', num2str(sfmin));
end

if ~subset
    colindx = 1:sf/sfmin:ncols;
    maptop = maplatlim(2);
    mapleft = maplonlim(1);
else
    colindx = collim(1):sf/sfmin:collim(2);
    maptop = corners(1,2) + dy*(rowlim(2)-1) + halfcell;
    mapleft = corners(1,1) + dx*(collim(1)-1) - halfcell;
end

scols = startprofiles(colindx);
cols = length(scols);

% Read from left to right of map
for n=1:cols
    fseek(fid,scols(n),'bof');
    fscanf(fid,'%s',[1 2]);						  % data element 1  (skip)
    profilenumrowscols = fscanf(fid,'%d',[1 2]);  % data element 2
    nrows = profilenumrowscols(1);
    switch ~subset
        case 1
            rowindx = 1:sf:nrows;
        otherwise
            rowindx = rowlim(1):sf:rowlim(2);
    end
    fscanf(fid,'%s',5);							  % data elements 3-5
    profile = fscanf(fid,'%d',1201);			  % data element 6
    if length(profile) ~= 1201
        fclose(fid);
        error(['map:' mfilename ':unableToReadFile'], ...
            'Could not read a whole profile. Perhaps the DEM does not conform to standard.');
    end
    map(:,n) = profile(rowindx); %#ok<AGROW>
end
fclose(fid);

cellsize = 1/celldim;
maplegend = [cellsize maptop mapleft];

%  --- Read Record Type C (Statistical Data) ---
%  Not implemented yet.

if nargout == 3
    Astruc = explainA(Astruc);
end

%--------------------------------------------------------------------------

function a = Adescription
% Data Set Identification (A) record contents

a( 1).length = 40;    a( 1).name = 'Quadrangle name';
a( 2).length = 40;    a( 2).name = 'Textual Info';
a( 3).length = 55;    a( 3).name = 'Filler';
a( 4).length = 1;     a( 4).name = 'Process Code';
a( 5).length = 1;     a( 5).name = 'Filler2';
a( 6).length = 3;     a( 6).name = 'Sectional Indicator';
a( 7).length = 4;     a( 7).name = 'MC origin Code';
a( 8).length = 1;     a( 8).name = 'DEM level Code'; 						a(8).type = '%6g';
a( 9).length = 1;     a( 9).name = 'Elevation Pattern Code'; 				a(9).type = '%6g';
a(10).length = 1;     a(10).name = 'Planimetric Reference System Code'; 	a(10).type = '%6g';
a(11).length = 1;     a(11).name = 'Zone';      							a(11).type = '%6g';
a(12).length = 15;    a(12).name = 'Projection Parameters'; 				a(12).type = '%24D';
a(13).length = 1;     a(13).name = 'Horizontal Units';  					a(13).type = '%6g';
a(14).length = 1;     a(14).name = 'Elevation Units';  						a(14).type = '%6g';
a(15).length = 1;     a(15).name = 'N sides To Bounding Box'; 				a(15).type = '%6g';
a(16).length = 8;     a(16).name = 'Bounding Box';     						a(16).type = '%24D';
a(17).length = 2;     a(17).name = 'Min Max Elevations';        			a(17).type = '%24D';
a(18).length = 1;     a(18).name = 'Rotation Angle';     					a(18).type = '%24D';
a(19).length = 1;     a(19).name = 'Accuracy Code';       					a(19).type = '%6g';
a(20).length = 3;     a(20).name = 'XYZ resolutions ';        				a(20).type = '%12E';
a(21).length = 2;     a(21).name = 'Nrows Cols';      						a(21).type = '%6g';

% Old format stops here

a(22).length = 1;     a(22).name = 'MaxPcontourInt';      					a(22).type = '%5g';
a(23).length = 1;     a(23).name = 'SourceMaxCintUnits';     				a(23).type = '%1g';
a(24).length = 1;     a(24).name = 'Smallest Primary';        				a(24).type = '%5g';
a(25).length = 1;     a(25).name = 'SourceMinCintUnits';      				a(25).type = '%1g';
a(26).length = 1;     a(26).name = 'Data Source Date';        				a(26).type = '%4g';
a(27).length = 1;     a(27).name = 'DataInspRevDate';       				a(27).type = '%4g';
a(28).length = 1;     a(28).name = 'InspRev Flag';
a(29).length = 1;     a(29).name = 'DataValidationFlag'; 	      			a(29).type = '%1g';
a(30).length = 1;     a(30).name = 'SuspectVoidFlag';    	    			a(30).type = '%2g';
a(31).length = 1;     a(31).name = 'Vertical Datum';       					a(31).type = '%2g';
a(32).length = 1;     a(32).name = 'Horizontal Datum';        				a(32).type = '%2g';
a(33).length = 1;     a(33).name = 'DataEdition';       					a(33).type = '%4g';
a(34).length = 1;     a(34).name = 'Percent Void';      					a(34).type = '%4g';

for i=1:length(a);
    if isempty(a(i).type)
        a(i).type = 'char';
    end
end

%--------------------------------------------------------------------------

function A = explainA(A)
% Fill in text associated with numeric codes in the USGS DEM A record

if ~isempty(A.ProcessCode)
    switch A.ProcessCode
        case '1' , A.ProcessCode = 'GPM';
        case '2' , A.ProcessCode = 'Manual Profile';
        case '3' , A.ProcessCode = 'DLG2DEM';
        case '4' , A.ProcessCode = 'DCASS';
    end
end

if ~isempty(A.ElevationPatternCode)
    switch A.ElevationPatternCode
        case 1, A.ElevationPatternCode = 'regular';
        case 2, A.ElevationPatternCode = 'random';
    end
end

if ~isempty(A.PlanimetricReferenceSystemCode)
    switch A.PlanimetricReferenceSystemCode
        case 0, A.PlanimetricReferenceSystemCode = 'Geographic';
        case 1, A.PlanimetricReferenceSystemCode = 'UTM';
        case 2, A.PlanimetricReferenceSystemCode = 'State Plane';
    end
end

if ~isempty(A.HorizontalUnits)
    switch A.HorizontalUnits
        case 0, A.HorizontalUnits = 'radians';
        case 1, A.HorizontalUnits = 'feet';
        case 2, A.HorizontalUnits = 'meters';
        case 3, A.HorizontalUnits = 'arc-seconds';
    end
end

if ~isempty(A.ElevationUnits)
    switch A.ElevationUnits
        case 1, A.ElevationUnits = 'feet';
        case 2, A.ElevationUnits = 'meters';
    end
end

if ~isempty(A.AccuracyCode)
    switch A.AccuracyCode
        case 0, A.AccuracyCode = 'unknown accuracy';
        case 1, A.AccuracyCode = 'accuracy information in record C';
    end
end

if ~isempty(A.SourceMaxCintUnits)
    switch A.SourceMaxCintUnits
        case 0, A.SourceMaxCintUnits = 'N.A.';
        case 1, A.SourceMaxCintUnits = 'feet';
        case 2, A.SourceMaxCintUnits = 'meters';
    end
end

if ~isempty(A.SourceMinCintUnits)
    switch A.SourceMinCintUnits
        case 1, A.SourceMinCintUnits = 'feet';
        case 2, A.SourceMinCintUnits = 'meters';
    end
end

if ~isempty(A.DataValidationFlag)
    switch A.DataValidationFlag
        case 0, A.DataValidationFlag = 'No validation performed';
        case 1, A.DataValidationFlag = 'TESDEM (record C added) no qualitative test (no DEM Edit System [DES] review)';
        case 2, A.DataValidationFlag = 'Water body edit and TESDEM run';
        case 3, A.DataValidationFlag = 'DES (includes water edit) no qualitiative test (no TESDEM)';
        case 4, A.DataValidationFlag = 'DES with record C added, qualitative and quantitative tests for level 1 DEM';
        case 5, A.DataValidationFlag = 'DES and TESTDEM qualitative and quantitative tests for levels 2 and 3 DEMs';
    end
end

if ~isempty(A.SuspectVoidFlag)
    switch A.SuspectVoidFlag
        case 0, A.SuspectVoidFlag = 'none';
        case 1, A.SuspectVoidFlag = 'suspect areas';
        case 2, A.SuspectVoidFlag = 'void areas';
        case 3, A.SuspectVoidFlag = 'suspect and void areas';
    end
end

if ~isempty(A.VerticalDatum)
    switch A.VerticalDatum
        case 1, A.VerticalDatum = 'local means sea level';
        case 2, A.VerticalDatum = 'National Geodetic Vertical Datum 1929 (NGVD 29)';
        case 3, A.VerticalDatum = 'North American Vertical Datum 1988 (NAVD 88)';
    end
end

if ~isempty(A.HorizontalDatum)
    switch A.HorizontalDatum
        case 1, A.HorizontalDatum = 'North American Datum 1927 (NAD27)';
        case 2, A.HorizontalDatum = 'World Geodetic System 1972 (WGS72)';
        case 3, A.HorizontalDatum = 'WGS84';
        case 4, A.HorizontalDatum = 'NAD83';
        case 5, A.HorizontalDatum = 'Old Hawaii Datum';
        case 6, A.HorizontalDatum = 'Puerto Rico Datum';
        case 7, A.HorizontalDatum = 'NAD 83 Provisional';
    end
end
