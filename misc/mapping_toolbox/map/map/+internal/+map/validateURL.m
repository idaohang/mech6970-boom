function validateURL(parameterName, serverURL)
%validateURL Validate serverURL input
%
%   validateURL(parameterName, serverURL) validates serverURL as valid URL
%   string.  parameterName is a string containing the name for the
%   serverURL parameter. serverURL is validated to be a row vector URL
%   string, containing the protocol, 'http, 'https', or 'file'.

% Copyright 2009 The MathWorks, Inc.
% $Revision: 1.1.6.1 $  $Date: 2009/05/14 17:05:43 $

% Validate serverURL.
validateattributes(serverURL, {'char'}, {'nonempty', 'vector', 'row'}, ...
    'validateURL', parameterName);

try
    url = java.net.URL(serverURL);
catch e
    error(e.identifier, ...
        'Expected %s to be a valid URL string.', parameterName);
end

% Validate the protocol. Note that in the error message, 'file' is not
% included since it is undocumented.
protocols = {'http', 'https','file'};
protocol = char(url.getProtocol);
assert(any(strcmpi(protocol, protocols)), ...
    'geoweb:validateURL:invalidProtocol', ...
    ['Expected %s to contain protocol type ', ...
    '''%s://'' or ''%s://''.'], ...
    parameterName, protocols{1}, protocols{2});

% Verify that a host has been provided.
protocol = [protocol '://'];
assert(numel(serverURL) > numel(protocol), ...
    'geoweb:validateURL:expectedHost', ...
    'Expected %s to contain a host name.', parameterName);
end
