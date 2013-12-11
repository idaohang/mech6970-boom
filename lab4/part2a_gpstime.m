%% Decode Data Bits
% This code is very shakey. This does do some detection on whether bit
% transitions happend exactly as they are supposed to, but if we have enough
% accel to throw off a PLL, then this is going to get really fucked up.
% 
% Right now it's best use is to exclude SV 20, leaving only 4 SV's

clear all; close all; clc
load part2a_tracking
fprintf('Decoding Data Bits with Thresholding\n')

data = zeros(acq.nsv,n_code_per);
data_IP_thold = 500;
data_trans_idx = cell(acq.nsv,1);
nav_msg = cell(acq.nsv,1);
ch_status = ones(acq.nsv,1);

word_len = 30;
subframe_len = 1500/5;
data_len = length(trackRes(1).IP);

for ch = 1:acq.nsv
  
  % thresholding
  data(ch,:) = trackRes(ch).IP > data_IP_thold;   %   --> 0/1
  % start tracking bits after the second uptick
  uptick_idx = find(diff(data(ch,:))==1);
  downtick_idx = find(diff(data(ch,:))==-1);
  % make sure that each uptick/downtick is spaced a multiple of 20 ms apart
  if ( ~all(rem(diff(uptick_idx),20)==0) ) || ( ~all(rem(diff(downtick_idx),20)==0) )
    fprintf(['problem on SV ' num2str(acq.svs(ch)) '. Skipping. \n'])
    ch_status(ch) = 0;
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



%% Get the index of the first subframe by finding preambles

fprintf('Finding preambles from IP data (Akos)\n')

TLM_starts = [];
search_start = 0;
while search_start+6000 < data_len
  TLM_starts(end+1,:) = findPreambles_(trackRes, acq, ch_status, search_start);
  search_start = search_start + 6000;
end

n_subframes = size(TLM_starts);
n_subframes = n_subframes(1);

%% Get the Z-Counts

TOW_bits = cell(size(TLM_starts));
crsr = cell(size(TLM_starts));
TLM_parity = zeros(size(TLM_starts));
TOW = zeros(size(TLM_starts));
HOW_parity_chk = zeros(size(TLM_starts));
TLM_parity_chk = zeros(size(TLM_starts));

% iterate over each SV
for ch = 1:acq.nsv
  if ~ch_status(ch), continue; end 
  for n = 1:n_subframes
    
    % get the actual bits
    TOW_bits{n,ch} = zeros(1,17);
    crsr{n,ch} = zeros(1,17);    
    TLM_parity(n,ch) = sign( trackRes(ch).IP( TLM_starts(n,ch) + 29*20 + 1 ) );
    for k = 0:16
      crsr{n,ch}(k+1) = TLM_starts(n,ch) + (word_len)*20 + k*20 + 1;
    end
    TOW_bits{n,ch} = sign(trackRes(ch).IP(crsr{n,ch}));
    TOW_bits{n,ch} = xor( TLM_parity(n,ch)>0 , TOW_bits{n,ch}>0 ) == 0;
    % this is the GPS TOW corresponding to TLM_starts
    TOW(n,ch) = (bin2dec(num2str(TOW_bits{n,ch}))-1)*6;
    
  end  
end


%% Plot IP and Data Bits

close all

for k = 1
  
  fh = figure;
%     subplot(4,1,1:2)
      plot(trackRes(k).IP); hold on
      for kk = 1:n_subframes
        plot([TLM_starts(kk,k) TLM_starts(kk,k)],[-5000 5000],'r','LineWidth',4)
      end
      grid on
      xlabel('Time (ms)');
      ylabel('IP');
      title(['In-phase Prompt for PRN' num2str(acq.svs(k))])
      
%     subplot(4,1,3)
%       plot(data(k,:),'LineWidth',4)
%       grid on; hold on
%       for kk = 1:length(data_trans_idx{k,1})
%         idx = data_trans_idx{k,1}(kk);
%         plot([idx idx],[-.5 1.5],'r')
%       end
%       title('Thresholded Data Bits')
%       xlabel('Time (ms)')
%       ylim([-1 2])
%     subplot(4,1,4)
%       stairs(nav_msg{k,1})
%       ylim([-1 2])
%       xlabel('Data Chip #')
%       title('Decoded Data Message (0/1)')
%       grid on

%   saveas(fh,['part2a_data_bits_prn_' num2str(acq.svs(k)) '.fig']);
  
end






%% End Matters

fprintf('saving data to `part2a_gpstime.mat`\n')
save part2a_gpstime


