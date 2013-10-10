function answer = isnumeral(s)
%ISNUMERAL True for numerals.
%   For a string S, ISNUMERAL(S) is true for numerals and false otherwise.
%
%   See also ISLETTER, ISSPACE.

%   Copyright 1996-2003 The MathWorks, Inc.
%   $Revision: 1.1.10.3 $ $Date: 2005/11/15 01:07:20 $

% This implementation works for ASCII, not sure about UNICODE.
answer = ('0' <= s) & (s <= '9');
