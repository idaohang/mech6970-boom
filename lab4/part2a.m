%% MECH 6970 Lab4, Part 2, (a) - Primary
% 
% Loads the results of the acquisition refinement
% 
% @author Robert Cofield
% 
tic
clear all; close all; clc
acq = load('part2a_narrow_ack.mat');

filename = ['..' filesep 'data' filesep 'GPS_Data_NordNav1e.sim'];
fileid = fopen(filename);

%% Constants

fL1 = 154*10.23e6; % L1 frequency,  1.5754e+09 Hz
fs = 16.3676e6; % sampling frequency
fIF = 4.1304e6; % intermediate frequency
fD = 50; % Data message chipping frequency, Hz

Tca_chip = 1.023e-6; % L1 C/A chip period
Tca = 1.0e-3; % L1 C/A sequence period
Td_chip = 1/fD;
d_frame_size = 1500; % 1500 bits per frame @ data msg chip rate
d_subframe_size = d_frame_size/5; % bits per subframe @ data msg chip rate
Td = d_frame_size*Td_chip; % data message frame period (1500 bits per frame), sec
Ts = 1/fs; % sampling period

preamble = bin2dec('10001011');  

% integration_period = Td/5; % seconds of data to read in - 1 subframe
% total_data_bytes = integration_period/Td_chip


%% Generate PRNs with shift found in Acquisition

% prns = genprn(acq.svs, 1023, [-1 1],
prns_aligned0 = zeros(nsv,N);
for sv = 1:acq.nsv
%   prns_aligned0(sv,:) = shift(prn
end

%% Decode Data Bits


preamble_found = false;
while ~preamble_found

  for svpi = 1:acq.nsv
    sv = present_svs(svpi);
    % get local variables for the present sv
    tau_ = acq.tau_samples(svpi);
    fdopp_ = acq.fdopp(svpi);
    feff_ = fdopp_ + fIF;

  end
end




%% End Matters
fclose(fileid);
clear fileid
save part2a
toc
