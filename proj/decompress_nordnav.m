clear all; close all; clc

data_dir = ['..' filesep 'data' filesep 'final_proj_data' filesep];

bin_filename = 'run_south.bin';
sim_filename = 'run_south.sim';
% bin_filename = 'run_north.bin';
% sim_filename = 'run_north.sim';

% bin_filename = 'run1_2min.bin';
% sim_filename = 'run1_2min.sim';

bin_fid = fopen([data_dir bin_filename],'r');
bin_info = dir([data_dir bin_filename])
sim_fid = fopen([data_dir sim_filename],'w');

out_stop_size = 40000; % bytes

[byte_in, nread] = fread(bin_fid, 1, 'uint8');
bytes_read = 0;
bytes_written = 0;

hwb = waitbar(0, ['Parsed 0 / ' num2str(bin_info.bytes) 'bytes']);

while nread % stop when nothing gets read in
  
  bytes_read = bytes_read+1;
  
  if ~rem(bytes_read,10000)
    waitbar(bytes_read/bin_info.bytes, hwb, ['Parsed ' num2str(bytes_read/1000) ' / ' num2str(bin_info.bytes/1000) ' kB']);
    if bytes_written >= out_stop_size
      break
    end
  end
  
  binstr = de2bi(byte_in);
  len = length(binstr);
  if len < 8 % pad with zeros
    binstr = [binstr zeros(1,8-len)];
  end
  
  for k = 1:2:8
    switch num2str(binstr(k:k+1))
      case '1  1'
        new_byte = -1;
      case '1  0'
        new_byte = -3;
      case '0  0'
        new_byte = 1;
      case '0  1'
        new_byte = 3;
      otherwise
        warning(['Uh oh @ read byte ' num2str(bytes_read)]);
        continue
    end
    fwrite(sim_fid,new_byte,'int8');
  end
  
  % get next byte
  [byte_in, nread] = fread(bin_fid, 1, 'uint8');
  
  bytes_written = bytes_written+1;
  
end

sim_info = dir([data_dir sim_filename])

close(hwb);
fclose('all');