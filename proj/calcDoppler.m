function dop=calcDoppler(user_v, sat_v, dir_v)
%user velocity vector, satelltie velocity vector, direction vector (user to
%satellite).  All in same reference frame.

c=299792458;

dir_v=dir_v./norm(dir_v);
v_rel=user_v-sat_v;

v_rel=(dot(v_rel,dir_v)./norm(dir_vel)).*(dir_v./norm(dir_v));

dop=v_rel/c;

%the change in frequency = dop * initial_frequency;
%OR
% new frequency=initial frequency*(1+dop)