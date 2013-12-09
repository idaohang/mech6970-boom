% open data file
filename = ['..' filesep 'data' filesep 'GPS_Data_NordNav1e.sim'];

% constants
freqL1 = 154*10.23e6;
sample_frequency = 16.3676e6;
intermediate_frequency = 4.1304e6;
integration_period = 1e-3;
gpsPi = 3.14159265359; 

% open file
fid = fopen(sprintf('%s',filename));

% number of bytes
bytes_to_read = round(sample_frequency*integration_period);

% read  1 millisecond chunk of data
signal1 = fread(fid,bytes_to_read,'int8')';

% read  another 1 millisecond chunk of data
signal2 = fread(fid,bytes_to_read,'int8')';

% close file when done
fclose(fid);


