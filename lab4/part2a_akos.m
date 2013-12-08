%% MECH 6970 Lab4, Part 2, (a) - Primary
% 
% Loads the results of the acquisition refinement
% 
% @author Robert Cofield
% 
tic
clear all; clc; 
acq = load('part2a_narrow_ack.mat');
akos_track = load('akos_trackingResults');

filename = ['..' filesep 'data' filesep 'GPS_Data_NordNav1e.sim'];
fileid = fopen(filename);

%% Constants

svidx_ = 1;

fL1 = 154*10.23e6; % L1 frequency,  1.5754e+09 Hz
fs = 16.3676e6; % sampling frequency
fIF = 4.1304e6; % intermediate frequency
fD = 50; % Data message chipping frequency, Hz
fCA = 1.023e6; % C/A code frequency, Hz

TL1 = 1/fL1; % L1 carrier period
Tca_chip = 1.023e-6; % L1 C/A chip period
Tca = 1.0e-3; % L1 C/A sequence period
Td_chip = 1/fD;
Ts = 1/fs; % sampling period
Tid = Tca; % integrate & dump period (sec) - use the C/A code period (1ms)
Tnav_chip = 20e-3; % period of a single nav msg chip
Tnav_frame = Tnav_chip*1500; % period of the entire nav msg frame
Tnav_subframe = Tnav_frame/5; % period of a nav msg subframe

stop_time = 1; % how long to go in sec

preamble = bin2dec('10001011');  

% DLL filter 2nd order
dll = struct(...
  'tau1', 0.0698469387755102,...
  'tau2', 0.37...
  );

% PLL filter 2nd Order
pll = struct(...
  'tau1', 0.000111755102040816,...
  'tau2', 0.0296...
  );

% % Look at PLL filter
% num = [4*pi*pll.z*pll.f , (2*pi*pll.f)^2];
% den = [1 , 4*pi*pll.z*pll.f ,  (2*pi*pll.f)^2];
% rlocus(tf(num,den,Tid))
% pause


%% Tracking

n_code_per = stop_time/Tid; % how many data points we'll have in the end

for ch = 1:acq.nsv
  channel(ch).PRN = acq.svs(ch);
  channel(ch).acquiredFreq = acq.fdopp(ch) + fIF;
  channel(ch).codePhase = round( 1023*fs/fCA - acq.tau_samples(ch));
  channel(ch).status = 'T';
end

settings = akos_track.settings;
settings.msToProcess = n_code_per;
settings.numberOfChannels = acq.nsv;

[trackingResults,trackingChannel] = tracking(fileid, channel, settings);

%% Decode Data Bits

data = zeros(acq.nsv,n_code_per);
data_IP_thold = 1000;
data_trans_idx = cell(acq.nsv,1);
nav_msg = cell(acq.nsv,1);

for ch = 1:acq.nsv
  % thresholding
  data(ch,:) = trackingResults(ch).I_P > data_IP_thold;  
  % start tracking bits after the second uptick
  uptick_idx = find(diff(data(ch,:))==1);
  downtick_idx = find(diff(data(ch,:))==-1);
  % make sure that each uptick is spaced a multiple of 20 ms apart
  if ( ~all(rem(diff(uptick_idx),20)==0) ) || ( ~all(rem(diff(downtick_idx),20)==0) )
    fprintf(['problem on SV ' num2str(acq.svs(ch)) '. Skipping. \n'])
    continue;
  end
  begin_up = uptick_idx(2)+1;
  begin_down = downtick_idx(2)+1;
  begin = min(begin_up,begin_down);
  data_trans_idx{ch,1}(1) = begin;
  cnt = 1; % how many data bits we have;
  while true
    idx = data_trans_idx{ch,1}(end);
    % save the value of the current bit
    nav_msg{ch,1}(end+1) = data(ch,idx);
    
    if data_trans_idx{ch,1}(end) + 20 > n_code_per % reached the end of the data stream
      break;
    end
    % save the index of the next data bit within the signal stream
    data_trans_idx{ch,1}(end+1) = idx+20; % save the next index
  end
  
end



%% Plot IP

close all

for k = 1:2
  
figure;
  subplot(3,1,1)
    plot(trackingResults(k).I_P)
    grid on
    xlabel('Time (ms)');
    ylabel('IP');
    title(['In-phase Prompt for PRN' num2str(acq.svs(k))])
  subplot(3,1,2)
    plot(data(k,:),'LineWidth',4)
    grid on; hold on
    for kk = 1:length(data_trans_idx{k,1})
      idx = data_trans_idx{k,1}(kk);
      plot([idx idx],[-.5 1.5],'r')
    end
    title('Thresholded Data Bits')
    xlabel('Time (ms)')
    ylim([-1 2])
  subplot(3,1,3)
    stairs(nav_msg{k,1})
    ylim([-1 2])
    xlabel('Data Chip #')
end




%% End Matters

fclose(fileid);
clear fileid
save part2a
toc
