%% MECH 6970 Lab4, Part 2, (a) - Primary
% 
% Loads the results of the acquisition refinement
% 
% @author Robert Cofield
% 
tic
clear all; clc; 
acq = load('part2a_narrow_ack.mat');

filename = ['..' filesep 'data' filesep 'GPS_Data_NordNav1e.sim'];
fileid = fopen(filename);

trackRes = akos_tracking(fileid, acq);


close(hwb);
fprintf('Tracking finished. Saving workspace as `part2aTrackRes.mat`\n')
save part2a_tracking


%% End Matters

fclose(fileid);
clear fileid
% save part2a
toc
