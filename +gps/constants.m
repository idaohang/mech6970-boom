%% GPS Constants

% Speed of light, (m/s)
c = 299792458;

% PRN bit positions to use in the G2 register for each PRN #
tap = [2 6;
       3 7;
       4 8;
       5 9;
       1 9;
       2 10;
       1 8;
       2 9;
       3 10;
       2 3;
       3 4;
       5 6;
       6 7;
       7 8;
       8 9;
       9 10;
       1 4;
       2 5;
       3 6;
       4 7;
       5 8;
       6 9;
       1 3;
       4 6;
       5 7;
       6 8;
       7 9;
       8 10;
       1 6;
       2 7;
       3 8;
       4 9;
       5 10;
       4 10;
       1 7;
       2 8;
       4 10];

% Signal frequencies
fCA = 1.023e6; % C/A code frequency, Hz
fPY = 10.23e6; % P(Y) code frequency, Ha
fL1 = 1575.42e6; % L1 frequency, Hz
fL2 = 1227.60e6; % L2 frequency, Hz

% Minimum received signal power
PC1 = -158.5; % L1 C/A code, dBW
PY1 = -161.5; % L1 P(Y) code, dBW
PY2 = -164.5; % L2 P(Y) code, dBW