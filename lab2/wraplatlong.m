function [lat,long]=wraplatlong(longin, latin)

latin=wrapTo360(latin);
if latin>180
    latin=latin-360;
end

if latin>90
    latin=90-abs(latin-90);
    longin=wrapTo360(longin+180);
elseif latin<-90
    latin=-90+abs(latin+90);
    longin=wrapTo360(longin-180);
end

lat=latin;
long=wrapTo360(longin+360);