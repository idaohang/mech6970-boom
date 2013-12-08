%% MECH 6970 Lab4 Part 1, b)
% Andrew Collins
% 
% Assumes Data file is in same folder, check name
% 
% NOTE !!! :
% The bit-level operations and array sizing for this 1ms analysis should not be
% exactly replicated for the 10ms part.
%   - in particular, the calculation of the upsample rate needs to be redone.
% 
clear all; close all; clc;
fprintf('Part I b)\n')
tic 
% try matlabpool 3; catch e, disp(e); end

%% constants

filename = ['..' filesep 'data' filesep 'GPS_Data_NordNav1e.sim'];

fL1 = 154*10.23e6; % L1 frequency,  1.5754e+09 Hz
fs = 16.3676e6; % sampling frequency
fIF = 4.1304e6; % intermediate frequency

Ts = 1/fs; % sampling period
integration_period = 1.0e-2; % Integration period changes to 10 ms
Tca = 1.023e-6; % L1 C/A code period
tau_chip_size = 0.5; % how many chips to use for tau

nfdopp = 25; % number of fdopp bins - 1
fdopp_bound = 5e3; % boundaries on either side of fIF to search for fdopp
PC1 = -158.5; % Power of L1 C/A code, dBW

fprintf('Integration Period: %f sec\n',integration_period)
fprintf('# of Bins in doppler freq search: %i\n',nfdopp+1)
fprintf('Search space for doppler freq: \n\t%f Hz  ---  %f Hz\n',-fdopp_bound,fdopp_bound)
fprintf('Tau: %f chips\n',tau_chip_size)


%% Read Nordnav Data

fid = fopen(filename);
bytes_to_read = round(fs*integration_period); % number of bytes
signal1 = fread(fid,bytes_to_read,'int8')'; % read  1 millisecond chunk of data
signal2 = fread(fid,bytes_to_read,'int8')'; % read another 1 millisecond chunk of data
fclose(fid); % close file when done

N = length(signal1);
noise1=6*randn(1,N); 
noise2=12*randn(1,N); 
signal1=signal1+noise1+noise2; 
T = 0:Ts:Ts*(N-1); % time from start corresponding to each epoch
upsample = N/1023; 
fprintf('Upsampling Ratio: %f\n',upsample)

clear filename fid  ans


%% Acquisition

% PRN replica
% generate PRNS for each satellite at the appropriate frequency
prn = cacode(4,upsample/10,10);


% Generate doppler stuff
dfdopp = 2*fdopp_bound/nfdopp; % delta fdopp, Hz
fprintf('Size of Bins in doppler freq search: %f Hz\n',dfdopp)
fdopp = linspace(-fdopp_bound, fdopp_bound, nfdopp+1); % fdopp in Hz
feff = fIF + fdopp; % effective frequency, Hz

tau_idx = tau_chip_size*round(upsample); % how many data indices correspond to tau
tau = [tau_idx:tau_idx:round(N),163676]; % vector of tau indices to use

y = cell(1,32); % each cell index contains an array for each sv
fdopp_soln = zeros(1,32);
fdopp_soln_idx = zeros(1,32);
tau_soln = zeros(1,32);
tau_soln_idx = zeros(1,32);
% parfor sv = 1:32 % generate correlation grid for each SV
for sv = 4
  y{1,sv} = zeros(length(tau),length(fdopp)); % signal replica for correlation
  for t_ = 1:length(tau) % loop over time shift values
    prn_shifted = shift(prn(1,:),tau(t_));
    for fd_ = 1:length(fdopp); % loop over doppler frequency values
      sin_ = imag(exp(1j*2*pi*feff(fd_)*T));
      cos_ = real(exp(1j*2*pi*feff(fd_)*T));
      I = signal1.*prn_shifted.*sin_;
      Q = signal1.*prn_shifted.*cos_;
      y{1,sv}(t_,fd_) = sum(I)^2 + sum(Q)^2;
    end
  end
  % find our answer
  [max_, fdopp_max_idx] = max(y{1,sv},[],2);
  [~, tau_max_idx] = max(max_);
  fdopp_soln_idx(sv) = fdopp_max_idx(tau_max_idx);
  fdopp_soln(sv) = fdopp(fdopp_soln_idx(sv));
  tau_soln_idx(sv) = tau_max_idx;
  tau_soln(sv) = tau(tau_soln_idx(sv));
end

fprintf('\nSOLUTION:\n')
fprintf('\tDoppler Shift: %f Hz\n',fdopp_soln(4))
fprintf('\tTau (Arrival Time): \n\t\t%f samples\n\t\t%f chips\n',tau_soln(4),tau_soln(4)/upsample)

%% Plot SV 4
    
fh = figure;
  surf(fdopp,tau,y{4})
  xlabel('Doppler Frequency (Hz)'); ylabel('Sample Shifts'); zlabel('Correlation');
  title('AutoCorrelation of Signal with Replica (10 ms), PRN #4');
saveas(fh, 'part1b.fig');


%% End matters
% try matlabpool close; catch e, disp(e); end
toc

