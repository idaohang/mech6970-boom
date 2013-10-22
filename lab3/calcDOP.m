function DOP = calcDOP(H)

% Calculate Dilutions of precision, given a Geometry Matrix
% 
% INPUTS:
%   H: Geometry Matrix
% OUTPUTS:
%   DOP: all the various dilution of precision information


H_ = inv(H'*H);

HDOP = sqrt(H_(1,1) + H_(2,2));
VDOP = sqrt(H_(3,3));
PDOP = sqrt(H_(1,1) + H_(2,2) + H_(3,3));
TDOP = sqrt(H_(4,4));
GDOP = sqrt(trace(H_));

DOP = [HDOP, VDOP, PDOP, TDOP, GDOP];

end