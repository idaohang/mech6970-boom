%% MECH 6970 Lab4 Part 1, a)
% Robert Cofield
% 
% 
clear all; close all; clc;
tic 
try matlabpool 2; catch e, disp(e); end

%% constants

filename = ['..' filesep 'data' filesep 'GPS_Data_NordNav1e.sim'];

fL1 = 154*10.23e6; % L1 frequency,  1.5754e+09 Hz
fs = 16.3676e6; % sampling frequency
fIF = 4.1304e6; % intermediate frequency

Ts = 1/fs; % sampling period
integration_period = 1.0e-3; % grab 1ms of data
Tca = 1.023e-6; % L1 C/A code period

nfdopp = 50; % number of fdopp bins - 1
fdopp_bound = 5e3; % boundaries on either side of fIF to search for fdopp
PC1 = -158.5; % Power of L1 C/A code, dBW


%% Read Nordnav Data

fid = fopen(filename);
bytes_to_read = round(fs*integration_period); % number of bytes
signal1 = fread(fid,bytes_to_read,'int8')'; % read  1 millisecond chunk of data
signal2 = fread(fid,bytes_to_read,'int8')'; % read another 1 millisecond chunk of data
fclose(fid); % close file when done

N = length(signal1);
T = 0:Ts:Ts*(N-1); % time from start corresponding to each epoch
upsample = N/1023; 

clear filename fid bytes_to_read ans


%% 

% PRN replica
% generate PRNS for each satellite at the appropriate frequency
prn = genprn(1:32, 1023, [-1 1], upsample);

% Generate doppler stuff
dfdopp = 2*fdopp_bound/nfdopp; % delta fdopp, Hz
fdopp = linspace(-fdopp_bound, fdopp_bound, nfdopp+1); % fdopp in Hz
feff = fIF - fdopp; % effective frequency, Hz

tau_idx = 0.5*upsample; % how many data indices correspond to tau
tau = tau_idx:tau_idx:1023*upsample; % vector of tau indices to use

y = cell(1,32);
parfor sv = 1:32 % generate correlation grid for each SV
  y{sv} = zeros(length(tau),length(fdopp)); % signal replica for correlation
  for t_ = 1:length(tau) % loop over time shift values
    for fd_ = 1:length(fdopp); % loop over doppler frequency values
      tau_ = tau(t_);
      fdopp_ = fdopp(fd_);
      f_eff = fIF + fdopp_; % effective frequency
      for n = 1:length(T); % summation over all time values
        hf_term = exp(1j*4*pi*f_eff*T(n));
        prn_shifted = shift(prn(32,:),n-tau_);
        prn_value = prn_shifted(1);
        y{sv}(t_,fd_) = y{sv}(t_,fd_) + ... %% !!! FINISH ME
          sqrt(2*PC1)*prn_value;
      end
      y{sv}(t_,fd_) = y{sv}(t_,fd_)/N; % average by dividing by number of epochs
    end
  end
end
    
  



%% End matters
try matlabpool close; catch e, disp(e); end
toc

