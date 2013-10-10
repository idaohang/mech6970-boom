function angmat = angledim(angmat,from,to)
%ANGLEDIM Convert angle units
%
%   ANGLEDIM has been replaced by four more specific functions:
%
%        fromRadians
%        fromDegrees
%        toRadians
%        toDegrees
%
%   but ANGLEDIM will be maintained for backward compatibility.  The
%   functions DEGTORAD, RADTODEG, and UNITSRATIO provide additional
%   alternatives.
%
%   angleOut = ANGLEDIM(angleIn, FROM, TO) converts angleIn from angle
%   units FROM to angle units TO.  FROM and TO may be either 'degrees'
%   or 'radians'.  They are case-insensitive and may be abbreviated.
%   angleIn and angleOut are arrays of class double, and size(angleOut)
%   matches size(angleIn).
%
%   Alternatives to ANGLEDIM
%   ------------------------
%   Because it must resolve both the input and output units, ANGLEDIM is
%   excessive for most applications.  In addition, it works only for
%   class double and it quietly discards the imaginary part of any
%   complex input.  Consider one of the following alternatives for
%   improved efficiency and generality:
%
%   If you are working from the command line, you can often replace
%   ANGLEDIM with DEGTORAD or RADTODEG.
%
%   If you are converting angle units within a script or function and you
%   know the values of both FROM and TO at the time of coding, then you can
%   also replace ANGLEDIM with DEGTORAD or RADTODEG.
%
%   If you know either FROM or TO at the time of coding, then you can
%   use fromRadians, fromDegrees, toRadians, or toDegrees.  Apply one of
%   the following transformations to your code:
%
%     angledim(angleIn,'radians',TO) --> fromRadians(TO,angleIn)
%
%     angledim(angleIn,'degrees',TO) --> fromDegrees(TO,angleIn)
%
%     angledim(angleIn,FROM,'radians') --> toRadians(FROM,angleIn)
%
%     angledim(angleIn,FROM,'degrees') --> toDegrees(FROM,angleIn)
%
%   Also note that the functions in the fromRadians family can convert
%   multiple variables in a single function call.  For example, you can
%   replace:
%
%     angle1 = angledim(angle1InRadians,'radians',TO);
%     angle2 = angledim(angle2InRadians,'radians',TO);
%
%   with:
%
%     [angle1,angle2] = fromRadians(TO,angle1InRadians,angle2InRadians);
%
%   If you do not know either FROM or TO at the time of coding, then you
%   can call UNITSRATIO to obtain the correct conversion factor, then
%   multiply the values of one or more variables.  For example, you can
%   replace:
%
%     angle1Out = angledim(angle1In, FROM, TO);
%     angle2Out = angledim(angle2In, FROM, TO);
%
%   with:
%
%     r = unitsratio(TO, FROM);
%     angle1Out = r * angle1In;
%     angle2Out = r * angle2In;
%
%   See also DEGTORAD, fromDegrees, fromRadians, RADTODEG, toDegrees,
%            toRadians, UNITSRATIO.

% Copyright 1996-2009 The MathWorks, Inc.
% $Revision: 1.12.4.8 $  $Date: 2009/03/30 23:38:05 $

if ~isa(angmat, 'double')
    id = sprintf('%s:%s:nonDoubleInput', getcomp, mfilename);
    error(id,'ANGMAT must be double.');
end

from = checkAngleStrings(from);
to   = checkAngleStrings(to);

% Convert complex input.
if ~isreal(angmat)
    angmat = real(angmat);
end

% Optimize slightly for no-op situation.
if strcmp(from,to)
    return
end

%  Find the appropriate string matches and transform the angles
switch from
    case 'degrees'
	      switch to
				case 'dm',        angmat = deg2dm(angmat);
				case 'dms',       angmat = deg2dms(angmat);
		        case 'radians',   angmat = angmat*pi/180;
	     end

	case 'radians'
	      switch to
		        case 'degrees',   angmat = angmat*180/pi;
				case 'dm',        angmat = rad2dm(angmat);
				case 'dms',       angmat = rad2dms(angmat);
	     end

	case 'dm'  
        % This case is stage 2 obsolete for R2007a. A warning will be
        % thrown in DMS2DEG, DM2DMS or DMS2RAD.
	      switch to
		        case 'degrees',   angmat = dms2deg(angmat);
				case 'dms',       dm2dms % Just throw error
				case 'radians',   angmat = dms2rad(angmat);
	     end

	case 'dms' 
        % This case is stage 2 obsolete for R2007a. A warning will be
        % thrown in DMS2DEG, DMS2DM or DMS2RAD.
	      switch to
		        case 'degrees',   angmat = dms2deg(angmat);
				case 'dm',        angmat = dms2dm(angmat);
				case 'radians',   angmat = dms2rad(angmat);
	     end
end

%--------------------------------------------------------------------------

function dm2dms

eid = sprintf('%s:%s:obsoleteDMS',getcomp,mfilename);
error(eid,...
    '''dm'' and ''dms'' angle encodings are obsolete.')

%--------------------------------------------------------------------------

function out = checkAngleStrings(in)

if ~ischar(in)
    eid = sprintf('%s:%s:nonStrInput', getcomp, mfilename);
    error(eid,'FROM and TO must be strings.');
end

in = lower(in);
switch in
  case {'radians' 'r' 'ra' 'rad' 'radi' 'radia' 'radian'}
    out = 'radians';
    
  case {'degrees' 'de' 'deg' 'degr' 'degre' 'degree'}
    out = 'degrees';
    
  case {'dm'}
    out = 'dm';
    
  case {'dms'}
    out = 'dms';
    
  otherwise
    eid = sprintf('%s:%s:invalidAngleUnits', getcomp, mfilename);
    error(eid,'FROM and TO must be either ''degrees'' or ''radians''.');
end
