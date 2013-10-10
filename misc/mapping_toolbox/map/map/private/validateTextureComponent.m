function [dataArgs, refmat, imageIndex, rules] = ...
    validateTextureComponent(mapfilename, dataArgs, displayType, rules)
%VALIDATETEXTURECOMPONENT Validate data for surface texture map
%
%   [dataArgs, R, imageIndex, RULES] = validateTextureComponent(
%   MAPFILNAME, dataArgs, displayType, RULES) validates the data in the
%   cell array dataArgs for use in surface as a texture mapped object.

% Copyright 2010 The MathWorks, Inc.
% $Revision: 1.1.6.3 $  $Date: 2010/11/17 11:24:58 $

checkmapnargin(2,4,rules.numArgs,mapfilename,displayType);
switch rules.numArgs
    case 2
        imageIndex = 1;
        if strcmp(displayType,'image')
            % (I,  R)
            % (BW, R)
            % (RGB,R)
            cmap = [];
            dataArgs{imageIndex} = checkImage( ...
                mapfilename, dataArgs{1}, cmap, 1, 2);           
        else
            % (Z,R)
            validateattributes(dataArgs{imageIndex}, ...
                {'numeric', 'logical'}, ...
                {'real','nonempty'}, mapfilename,'Z', 1);
            sz = size(dataArgs{imageIndex}, 3);
            assert(sz == 3 || sz == 1, ...
                'map:validateTextureComponent:expected2Dor3D', ...
                ['Function %s expected its first input, Z, to be ', ...
                'size M-by-N or M-by-N-by-3.'], upper(mapfilename));
            if islogical(dataArgs{imageIndex})
                dataArgs{imageIndex} = double(dataArgs{imageIndex});
            end
        end
        R_pos = 2;
        R = dataArgs{R_pos};
        rasterSize = size(dataArgs{imageIndex});
        refmat = checkRefObj(mapfilename, R, rasterSize, R_pos);
        
    case 3
        refmat = [];
        R = [];
        if strcmp(displayType,'image')
            if rules.isImageC1C2
                % (C1, C2, I)
                % (C1, C2, BW)
                % (C1, C2, RGB)
                cmap = [];
                imageIndex = 3;
                dataArgs{imageIndex} = checkImage( ...
                    mapfilename, dataArgs{imageIndex}, cmap, 3, 2);
                require2D = false;
                checkGeolocatedDataGrid(dataArgs{1}, dataArgs{2}, ...
                   dataArgs{imageIndex}, rules, mapfilename, require2D);
            else
                % (X, CMAP, R) with mapped axes
                cmap = dataArgs{2};
                imageIndex = 1;
                dataArgs{imageIndex} = checkImage( ...
                    mapfilename, dataArgs{imageIndex}, cmap, 1, 2);
                R_pos = 3;
                R = dataArgs{R_pos};
                rasterSize = size(dataArgs{imageIndex});
                refmat = checkRefObj(mapfilename, R, rasterSize, R_pos);
            end
        else
            % (C1, C2, Z)
            imageIndex = 3;
            rules.isImageC1C2 = true;
            require2D = false;
            checkGeolocatedDataGrid(dataArgs{1}, dataArgs{2}, ...
                 dataArgs{imageIndex}, rules, mapfilename, require2D);
        end
        
    case 4
        % (C1, C2, X, CMAP)
        refmat = [];
        R = [];
        cmap = dataArgs{4};
        imageIndex = 3;
        dataArgs{imageIndex} = checkImage( ...
            mapfilename, dataArgs{imageIndex}, cmap, 3, 4);
        require2D = false;
        checkGeolocatedDataGrid(dataArgs{1}, dataArgs{2}, ...
            dataArgs{imageIndex}, rules, mapfilename, require2D);
end

% If R is a raster referencing object and the RasterInterpretation is
% 'postings', then issue an error.
if isobject(R) && strcmp(R.RasterInterpretation, 'postings')
    error('map:rastershow:postingsWithTexture', ...
        ['The raster referencing object, %s, is defined with raster ', ...
        'interpretation set to ''%s'' which is not supported when the ''%s'' ' ...
        'parameter is set to ''%s''. To display this data, ', ...
        'set the ''%s'' parameter to ''%s'', ''%s'', or ''%s''.'], ...
        'R', 'postings', 'DisplayType', displayType, 'DisplayType', ...
        'surface', 'mesh', 'contour');
end
    
% Convert to double if required.
isDoubleOrUint8 = isa(dataArgs{imageIndex}, 'double') ...
    || isa(dataArgs{imageIndex}, 'uint8');
if ~isDoubleOrUint8
    dataArgs{imageIndex} = convertToDouble(dataArgs{imageIndex});
end

%--------------------------------------------------------------------------

function img = convertToDouble(img)
% convertToDouble takes an image as input, and returns an array of class
% double.  If the input image is of class double, the output array is
% identical to it.  If the input image is not double, convertToDouble
% returns an array of class double, rescaling the data as necessary,
% according to the following table:
%
%     Data class (numeric or logical)         Operation
%     -------------------------------         ---------
%     double                                  none
%     uint8,uint16,uint32,uint64              rescaled, cast to double
%     single,int8,int16,int32,int64,logical   cast to double

if ~isa(img, 'double')
    if isa(img, 'uint8') || isa(img, 'uint16') ...
            || isa(img, 'uint32') || isa(img, 'uint64')
        classType = class(img);
        img = double(img) / double(intmax(classType));
    else
        img = double(img);
    end
end
