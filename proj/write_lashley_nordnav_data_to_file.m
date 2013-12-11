clear all; close all; clc
load /media/OneTouch' 4'/matlab' code'/NordNav_IF_data.mat
fid = fopen('../data/final_proj_data/lashley_nordnav_data.sim','w')
for k = 1:length(A)
  fwrite(fid,A(k),'int8');
end
fclose(fid);

fid = fopen('../data/final_proj_data/lashley_nordnav_data.sim','r')
B = fread(fid, length(A), 'int8');
all(A==B)
fclose(fid);