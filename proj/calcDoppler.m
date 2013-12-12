function dop=calcDoppler(user_pos, course, user_vel, sat_v, sat_pos)
%user position(lat,long,alt), LLA, course, user velocity magnitude (along course), satelltie position vector, direction vector (user to
%satellite).  All in same reference frame.
lat=user_pos(1);
long=user_pos(2);

c=299792458;

vel_E=cos(course)*user_vel;
vel_N=sin(course)*user_vel;
vel_U=0;
vel_ENU=[vel_E,vel_N,vel_U];

user_pos_ecef = wgslla2xyz(user_pos(1), user_pos(2), user_pos(3));
dir_v=sat_pos-user_pos_ecef';

user_v=rotenu2xyz(vel_ENU', lat, long);

dir_v=dir_v./norm(dir_v);
v_rel=user_v-sat_v';

v_rel=(dot(v_rel,dir_v)./norm(dir_v)).*(dir_v./norm(dir_v));

dop=norm(v_rel)/c;
freqL1 = 154*10.23e6;

dop=dop*freqL1;
%the change in frequency = dop * initial_frequency;
%OR
% new frequency=initial frequency*(1+dop)