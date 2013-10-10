function schema
%SCHEMA Definition for Text class

% Copyright 1996-2008 The MathWorks, Inc.
% $Revision: 1.1.6.5 $  $Date: 2008/10/26 14:25:55 $

pkg = findpackage('MapGraphics');
c = schema.class(pkg,'Text');

p = schema.prop(c,'hText','MATLAB array');
p.AccessFlags.PrivateGet = 'on';
p.AccessFlags.PrivateSet = 'on';
p.AccessFlags.PublicGet  = 'on';
p.AccessFlags.PublicSet  = 'off';
