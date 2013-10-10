function validateRasterReference(R, classes, func_name, var_name, arg_pos)
%validateRasterReference Validate referencing matrix or object
%
%   validateRasterReference(R, CLASSES, FUNC_NAME, VAR_NAME, ARG_POS)
%   ensures that R is either a valid referencing matrix or a scalar object
%   matching one of the class types specified in the cell array CLASSES.
%
%   Examples
%   --------
%   % Wrong size
%   clear R; R(4,5,6) = 0;
%   internal.map.validateRasterReference(R, {'spatialref.GeoRasterReference', ...
%       'spatialref.MapRasterReference'}, 'MY_FUNC', 'MY_VAR', 1)
%
%   % 3-by-2, but wrong class
%   R = zeros([3 2], 'uint8');
%   internal.map.validateRasterReference(R, ...
%       {'spatialref.MapRasterReference'}, 'MY_FUNC', 'MY_VAR', 1)
%
%   % Non-finite 3-by-2 double
%   R = [0 1; 1 0; -Inf 200];
%   internal.map.validateRasterReference(R, ...
%       {'spatialref.MapRasterReference'}, 'MY_FUNC', 'MY_VAR', 1)
%
%   % 3-by-2 complex double
%   R = [0 1; 1i 0; -100 200];
%   internal.map.validateRasterReference(R, ...
%       {'spatialref.MapRasterReference'}, 'MY_FUNC', 'MY_VAR', 1)
%
%   % Non-scalar object
%   R = georasterref(); R(2) = georasterref();
%   internal.map.validateRasterReference(R, ...
%       {'spatialref.GeoRasterReference'}, 'MY_FUNC', 'MY_VAR', 1)
%
%   % Unexpected class
%   R = georasterref();
%   internal.map.validateRasterReference(R, ...
%       {'spatialref.MapRasterReference'}, 'MY_FUNC', 'MY_VAR', 1)
%
%   % Both classes accepted: no exception thrown
%   R = maprasterref();
%   internal.map.validateRasterReference(R, {'spatialref.GeoRasterReference', ...
%       'spatialref.MapRasterReference'}, 'MY_FUNC', 'MY_VAR', 1)
%
%   % No classes are accepted but an object is input
%   % (Validate strictly as referencing matrix.)
%   R = georasterref();
%   internal.map.validateRasterReference(R, {}, 'MY_FUNC', 'MY_VAR', 1)
%
%   % No classes are accepted and input is a valid referencing matrix
%   R = [0 1; -1 0; 100 200];
%   internal.map.validateRasterReference(R, {}, 'MY_FUNC', 'MY_VAR', 1)

% Copyright 2010 The MathWorks, Inc.
% $Revision: 1.1.6.1 $  $Date: 2010/10/11 14:47:33 $

if ~isobject(R) || isempty(classes)
    % Validate referencing matrix. It must be a 3-by-2 matrix of real-valued
    % finite doubles.
    if ~isequal(size(R),[3,2])
        sizestr = sprintf('%dx', size(R));
        sizestr(end) = [];
        msg2 = sprintf('Instead its size was %s.', sizestr);
        if isempty(classes)
            throwAsCaller(MException('map:validate:expectedReferencingMatrix', ...
                'Function %s expected input number %d, %s, to a 3-by-2 referencing matrix. %s', ...
                func_name, arg_pos, var_name, msg2))
        else
            throwAsCaller(expectedReferencingMatrixOrObject( ...
                classes, func_name, var_name, arg_pos, msg2))
        end
    end
    
    % If we reach this line, we know for certain that R is 3-by-2. We'll
    % assume that a referencing matrix was intended (rather than an
    % object) and phrase the error messages accordingly.
    try
        validateattributes(R, {'double'} ,{'real','finite'}, ...
            func_name, var_name, arg_pos)
    catch exception
        mnemonic = extract_mnemonic(exception.identifier, func_name);
        switch mnemonic
            
            case 'invalidType'
                throwAsCaller(MException('map:validate:expectedClassDoubleRefmat', ...
                    'Function %s expected input number %d, %s, a 3-by-2 referencing matrix, to be class %s.', ...
                    func_name, arg_pos, var_name, 'double'))
                
            case 'expectedFinite'
                throwAsCaller(MException('map:validate:expectedFiniteRefmat', ...
                    'Function %s expected input number %d, %s, a 3-by-2 referencing matrix, to contain finite values.', ...
                    func_name, arg_pos, var_name))
                
            case 'expectedReal'
                throwAsCaller(MException('map:validate:expectedRealRefmat', ...
                    'Function %s expected input number %d, %s, a 3-by-2 referencing matrix, to contain real values.', ...
                    func_name, arg_pos, var_name))
                
            otherwise
                % We don't expect to reach this line.
                rethrow(exception)
        end
    end
else
    % Validate scalar object. Its type must match one of the classes
    % listed in CLASSES.
    try
        validateattributes(R, classes, {'scalar'}, func_name, var_name, arg_pos)
    catch exception
        mnemonic = extract_mnemonic(exception.identifier, func_name);
        switch mnemonic
            
            case 'expectedScalar'
                sizestr = sprintf('%dx', size(R));
                sizestr(end) = [];
                msg2 = sprintf('Instead its size was %s.', sizestr);
                throwAsCaller(expectedReferencingMatrixOrObject( ...
                    classes, func_name, var_name, arg_pos, msg2))
                
            case 'invalidType'
                msg2 = sprintf('Instead its type was: %s.', class(R));
                throwAsCaller(expectedReferencingMatrixOrObject( ...
                    classes, func_name, var_name, arg_pos, msg2))
                
            otherwise
                % We don't expect to reach this line.
                rethrow(exception)
        end
    end
end

%-----------------------------------------------------------------------------

function exception = expectedReferencingMatrixOrObject( ...
    classes, func_name, var_name, arg_pos, msg2)
% Construct MException describing expected input.

classlist = sprintf('%s, ', classes{:});
classlist(end-1:end) = [];

exception = MException('map:validate:expectedReferencingMatrixOrObject', ...
    'Function %s expected input number %d, %s, to be either a 3-by-2 referencing matrix or a scalar instance of one of these classes:\n\n  %s\n\n%s', ...
    func_name, arg_pos, var_name, classlist, msg2);

%-----------------------------------------------------------------------------

function mnemonic = extract_mnemonic(identifier, func_name)
n = numel(['MATLAB:' func_name ':']);
mnemonic = identifier(n+1:end);
