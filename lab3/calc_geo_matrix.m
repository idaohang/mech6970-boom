function G=calc_geo_matrix(sv_pos,user_pos)


[sv_num,n]=size(sv_pos);
for i=1:sv_num
unit_vect=(sv_pos(i,:)-user_pos)'/norm(sv_pos(i,:)-user_pos);
G(i,:)=[-unit_vect(1),-unit_vect(2),-unit_vect(3), 1];
end

