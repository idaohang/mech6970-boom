function schema
%SCHEMA Definition for Arrow class

% Copyright 1996-2008 The MathWorks, Inc.
% $Revision: 1.1.6.3 $  $Date: 2008/10/26 14:25:36 $

pkg = findpackage('MapGraphics');
c = schema.class(pkg,'ArrowHead');

p = schema.prop(c,'hPatch','MATLAB array');
p.AccessFlags.PrivateGet = 'on';
p.AccessFlags.PrivateSet = 'on';
p.AccessFlags.PublicGet  = 'on';
p.AccessFlags.PublicSet  = 'on';

p = schema.prop(c,'line','MATLAB array');
p.AccessFlags.PrivateGet = 'on';
p.AccessFlags.PrivateSet = 'on';
p.AccessFlags.PublicGet  = 'on';
p.AccessFlags.PublicSet  = 'on';

p = schema.prop(c,'ArrowWidth','double');
p.AccessFlags.PrivateGet = 'on';
p.AccessFlags.PrivateSet = 'on';
p.AccessFlags.PublicGet  = 'on';
p.AccessFlags.PublicSet  = 'on';
p.FactoryValue = 16/3;

p = schema.prop(c,'ArrowHeight','double');
p.AccessFlags.PrivateGet = 'on';
p.AccessFlags.PrivateSet = 'on';
p.AccessFlags.PublicGet  = 'on';
p.AccessFlags.PublicSet  = 'on';
p.FactoryValue = 24;
