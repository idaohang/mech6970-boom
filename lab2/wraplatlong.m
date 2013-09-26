function [lat,long]=wraplatlong(latin,longin)

if latin>90
    latin=(latin-90)+(-90);
    longin=wrapTo360(long+180);
elseif latin<-90
    
    
end
    lat=latin;
    long=longin;