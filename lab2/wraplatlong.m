function [lat,long]=wraplatlong(longin, latin)

if latin>90
    while latin>90
    latin=(latin-180);
    longin=wrapTo360(longin+180);
    end
elseif latin<-90
    while latin<-90
    latin=(latin+180);
    longin=wrapTo360(longin-180);
    end
end
    lat=latin;
    long=longin;