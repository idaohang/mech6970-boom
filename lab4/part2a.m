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

% PLL filter 2nd order
K_dll = 1.01;
a_dll = 0.5;


%% Decode Data Bits

len = stop_time/Tid; % how many data points we'll have in the end

nbyte_id = round(fs*Tid); % number of bytes in the integrate & dump period
upsample = nbyte_id/1023; % samples per CA code chip
EL_shift = 0.5*upsample; % number of chips to shift the code to get early & late
prns_orig = genprn(acq.svs, 1023, [-1 1], upsample);

carry_on = true;
nav_msg = zeros(acq.nsv,len); % bits from the nav msgs stored here
I = struct('E',zeros(acq.nsv,len),'P',zeros(acq.nsv,len),'L',zeros(acq.nsv,len));
Q = struct('E',zeros(acq.nsv,len),'P',zeros(acq.nsv,len),'L',zeros(acq.nsv,len));

signal_time = zeros(1,len);
signal_time_ = -Ts;
data_idx = 0;

% get initial estimates for the calculated values
tau_chips = acq.tau_samples; % code phase 
feff = acq.fdopp + fIF;
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
    prnP = shift(prns_orig(s,:), tau_chips(s)); % prompt code
    % get early and late PRNs
    prnE = shift(prnP, EL_shift);
    prnL = shift(prnP, -EL_shift);
    % calculate the carrier
    sin_ = imag(exp(1j*2*pi*feff(s)*signal_time_));
    cos_ = real(exp(1j*2*pi*feff(s)*signal_time_));
    
    % Integrate & Dump
    I.E(s,k) = sum(signal.*prnE.*sin_);
    I.P(s,k) = sum(signal.*prnP.*sin_);
    I.L(s,k) = sum(signal.*prnL.*sin_);
    Q.E(s,k) = sum(signal.*prnE.*cos_);
    Q.P(s,k) = sum(signal.*prnP.*cos_);
    Q.L(s,k) = sum(signal.*prnL.*cos_);
    
    % calculate the code error
    early = sqrt( I.E(s,k)^2 + Q.E(s,k)^2 );
    late  = sqrt( I.L(s,k)^2 + Q.L(s,k)^2 );
    prn_phase_err_ = (early-late)/(early+late);
    
    % calculate the phase error
    car_phase_err(s,k) = atan(Q.P(s,k)/I.P(s,k)) / (2*pi); % Hz
    
    % filter errors
    prn_phase_err(s,k+1) = (K_dll*a_dll+1)/(1-K_dll) * prn_phase_err_;
    
    
  end
  
  % save time of this cumulative step as the time at the end of this
  signal_time(k) = signal_time_(end);
  
end


%% Plotting

figh = figure;
  plot(signal_time,I.P(1,:)); grid on
  title(['Inphase Prompt - SV#' num2str(acq.svs(1))])
  xlabel('Signal Time (s)'); ylabel('IP');
  
figh = figure;
  plot(signal_time, prn_phase_err(1,1:end-1)); grid on
  title(['Filtered Code Phase Error - SV#' num2str(acq.svs(1))]);
  xlabel('Signal Time (s)'); ylabel('e');

figh = figure;
  plot(signal_time, car_phase_err(1,1:end)); grid on
  title(['Carrier Phase Error - SV#' num2str(acq.svs(1))])
  xlabel('Signal Time (s)'); ylabel('f');
  

%% End Matters

fclose(fileid);
clear fileid
save part2a
toc
