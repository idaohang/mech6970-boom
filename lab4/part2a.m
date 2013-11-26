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
Ts = 1/fs; % sampling period
Tid = Tca; % integrate & dump period (sec) - use the C/A code period
Tnav_chip = 20e-3; % period of a single nav msg chip
Tnav_frame = Tnav_chip*1500; % period of the entire nav msg frame
Tnav_subframe = Tnav_frame/5; % period of a nav msg subframe

stop_time = Tnav_subframe;

preamble = bin2dec('10001011');  

% DLL filter 2nd order
K_dll = 0.5;
a_dll = -0.5;

% PLL filter 2nd Order
K_pll = 0.5;
a_pll = -0.5;

% % Look at PLL filter
% num = K_pll*[1 a_pll]
% den = [1 -1]
% rlocus(tf(num,den,Tid))
% pause


%% Decode Data Bits

len = stop_time/Tid; % how many data points we'll have in the end

nbyte_id = round(fs*Tid); % number of bytes in the integrate & dump period
upsample = nbyte_id/1023; % samples per CA code chip

EL_shift = 0.5*upsample; % number of chips to shift the code to get early & late
prns_orig = genprn(acq.svs, 1023, [-1 1], upsample);

prns_ = zeros(size(prns_orig)); % stores current PRN sequence
% calculate initial shifted PRNs based on acquisition data
for s = 1:acq.nsv
  prns_(s,:) = shift(prns_orig(s,:),acq.tau_samples(s));
end

nav_msg = zeros(acq.nsv,len); % bits from the nav msgs stored here
I = struct('E',zeros(acq.nsv,len),'P',zeros(acq.nsv,len),'L',zeros(acq.nsv,len));
Q = struct('E',zeros(acq.nsv,len),'P',zeros(acq.nsv,len),'L',zeros(acq.nsv,len));

signal_time = zeros(1,len);
signal_time_ = -Ts;
data_idx = 0;

% get initial estimates for the calculated values
tau_chips = zeros(acq.nsv,len);
tau_chips(:,1) = acq.tau_samples; % code phase 
car_freq = zeros(acq.nsv,len);
car_freq(:,1) = acq.fdopp + fIF;
prn_phase_err = zeros(acq.nsv,len);
car_phase_err = zeros(acq.nsv,len);

for k = 1:len
  
  % the index that the end of the new integrated&dumped signal will have in the
  % total signal
  data_idx = data_idx+1;
  
  % load data of length equal to the I&D period
  signal = fread(fileid,nbyte_id,'int8')';
  % find the time (sec) corresponding to the current signal chunk
  signal_time_ = signal_time_(end)+Ts:Ts:signal_time_(end)+Ts*nbyte_id;
  
%   for svi = 1:acq.nsv % iterate over present SVs
  for s = 1 % just over SV 4

    % % Delay Lock Loop
    
    % shift the original PRN sequence by the current code phase
    prnP = prns_(s,:); % prompt code
    % get early and late PRNs
    prnE = shift(prnP, EL_shift);
    prnL = shift(prnP, -EL_shift);
    % calculate the carrier
    sin_ = imag(exp(1j*2*pi*car_freq(s,k)*signal_time_));
    cos_ = real(exp(1j*2*pi*car_freq(s,k)*signal_time_));
    
    % Integrate & Dump
    I.E(s,k) = sum(signal.*prnE.*sin_);
    I.P(s,k) = sum(signal.*prnP.*sin_);
    I.L(s,k) = sum(signal.*prnL.*sin_);
    Q.E(s,k) = sum(signal.*prnE.*cos_);
    Q.P(s,k) = sum(signal.*prnP.*cos_);
    Q.L(s,k) = sum(signal.*prnL.*cos_);
    
    early = sqrt( I.E(s,k)^2 + Q.E(s,k)^2 );
    late  = sqrt( I.L(s,k)^2 + Q.L(s,k)^2 );
    
    % calculate the code error
    prn_phase_err_ = (early-late)/(early+late);
    
    % calculate the phase error
    car_phase_err_ = atan(I.P(s,k)/Q.P(s,k)); % rad
    
    % Do discrete filtering if enough data already present
    if k>1
      % filter carrier phase error
      car_phase_err(s,k) = (K_pll*a_pll+1)/(1-K_pll)*car_phase_err(s,k-1) * car_phase_err_; % rad
      % filter code phase errors
      prn_phase_err(s,k) = (K_dll*a_dll+1)/(1-K_dll)*prn_phase_err(s,k-1) * prn_phase_err_; % chips ???
    else
      car_phase_err(s,k) = car_phase_err_; % rad
      prn_phase_err(s,k) = prn_phase_err_/Ts; % transfer code phase to chips from sec
    end    
    
    % save new values for code phase & carrier frequency
%     tau_chips(s,k) = prn_phase_err(s,k);
    prns_(s,:) = shift(prns_(s,:),round(prn_phase_err(s,k)));
    car_freq(s,k) = car_freq(s,k) + car_phase_err(s,k)/(2*pi); % convert car freq to Hz
    
    % propagate carrier frequency to next time step
    if k~=len, car_freq(s,k+1) = car_freq(s,k); end
    
  end
  
  % save time of this cumulative step as the time at the end of this
  signal_time(k) = signal_time_(end);
  
end


%% Plotting

% Plot the actual IP
figh = figure;
  plot(signal_time,I.P(1,1:len)); grid on
  title(['Inphase Prompt - SV#' num2str(acq.svs(1))])
  xlabel('Signal Time (s)'); ylabel('IP');
  
% % Take the FFT of IP
L = length(I.P(1,:));
NFFT = 2^nextpow2(L);
Y = fft(I.P(1,:),NFFT)/L;
f = fs/2*linspace(0,1,NFFT/2+1);
figh = figure;
  plot(f,2*abs(Y(1:NFFT/2+1))); grid on
  title('Single Sided Amplitude Spectrum of IP')
  xlabel('Frequency (Hz)'); ylabel('| Y(F) |')

% % Plot the Code Phase error
figh = figure;
  plot(signal_time, prn_phase_err(1,1:len)); grid on
  title(['Filtered Code Phase Error - SV#' num2str(acq.svs(1))]);
  xlabel('Signal Time (s)'); ylabel('e');

% Plot Carrier Phase Error and Frequency
figh = figure;
  subplot(2,1,1)
    plot(signal_time, car_phase_err(1,1:len)); grid on
    axh1 = gca;
    title(['Carrier Phase Error - SV#' num2str(acq.svs(1))])
    xlabel('Signal Time (s)'); ylabel('f');
  subplot(2,1,2)
    plot(signal_time, car_freq(1,1:len)); grid on
    axh2 = gca;
    title('Carrier Frequency'); ylabel('f_{car} (Hz)'); xlabel('Signal Time (s)');
  linkprop([axh1 axh2], {'XLim'});
    
%% End Matters

fclose(fileid);
clear fileid
save part2a
toc
