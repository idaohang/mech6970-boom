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


%% Plot IP and Data Bits

close all

for k = 1:2
  
  fh = figure;
    subplot(4,1,1:2)
      plot(trackRes(k).IP);
      grid on
      xlabel('Time (ms)');
      ylabel('IP');
      title(['In-phase Prompt for PRN' num2str(acq.svs(k))])
    subplot(4,1,3)
      plot(data(k,:),'LineWidth',4)
      grid on; hold on
      for kk = 1:length(data_trans_idx{k,1})
        idx = data_trans_idx{k,1}(kk);
        plot([idx idx],[-.5 1.5],'r')
      end
      title('Thresholded Data Bits')
      xlabel('Time (ms)')
      ylim([-1 2])
    subplot(4,1,4)
      stairs(nav_msg{k,1})
      ylim([-1 2])
      xlabel('Data Chip #')
      title('Decoded Data Message (0/1)')
      grid on

  saveas(fh,['part2a_data_bits_prn_' num2str(acq.svs(k)) '.fig']);
  
end

%% Get the index of the first subframe by finding preambles
fprintf('Finding preambles from IP data (Akos)\n')

% Preamble search can be delayed to a later point in the tracking results 
% to avoid noise due to tracking loop transients  
searchStartOffset = 0; 
 
%--- Initialize the firstSubFrame array ----------------------------------- 
firstSubFrame = zeros(1, acq.nsv); 
 
%--- Generate the preamble pattern ---------------------------------------- 
preamble_bits = [1 -1 -1 -1 1 -1 1 1]; 
 
% "Upsample" the preamble - make 20 vales per one bit. The preamble must be 
% found with precision of a sample. 
preamble_ms = kron(preamble_bits, ones(1, 20)); 
 
%--- Make a list of channels excluding not tracking channels -------------- 
activeChnList = find(ch_status);
 
%=== For all tracking channels ... 
for channelNr = activeChnList'
 
% Correlate tracking output with preamble ================================ 
    % Read output from tracking. It contains the navigation bits. The start 
    % of record is skiped here to avoid tracking loop transients. 
    %bits = trackResults(channelNr).I_P(1 + searchStartOffset : end); 
    bits = trackRes(channelNr).IP(1+searchStartOffset:end);
 
    % Now threshold the output and convert it to -1 and +1  
    bits(bits > 0)  =  1; 
    bits(bits <= 0) = -1; 
 
    % Correlate tracking output with the preamble 
    tlmXcorrResult = xcorr(bits, preamble_ms); 
 
% Find all starting points off all preamble like patterns ================ 
    clear index 
    clear index2 
 
    xcorrLength = (length(tlmXcorrResult) +  1) /2; 
 
    %--- Find at what index/ms the preambles start ------------------------ 
    index = find(... 
        abs(tlmXcorrResult(xcorrLength : xcorrLength * 2 - 1)) > 153)' + ... 
        searchStartOffset; 
 
% Analyze detected preamble like patterns ================================ 
    for i = 1:size(index) % For each occurrence 
 
        %--- Find distances in time between this occurrence and the rest of 
        %preambles like patterns. If the distance is 6000 milliseconds (one 
        %subframe), the do further verifications by validating the parities 
        %of two GPS words 
         
        index2 = index - index(i); 
 
        if (~isempty(find(index2 == 6000))) 
 
            %=== Re-read bit vales for preamble verification ============== 
            % Preamble occurrence is verified by checking the parity of 
            % the first two words in the subframe. Now it is assumed that 
            % bit boundaries a known. Therefore the bit values over 20ms are 
            % combined to increase receiver performance for noisy signals. 
            % in Total 62 bits mast be read : 
            % 2 bits from previous subframe are needed for parity checking; 
            % 60 bits for the first two 30bit words (TLM and HOW words). 
            % The index is pointing at the start of TLM word. 
            bits = trackRes(channelNr).IP(index(i)-40 : ... 
                                               index(i) + 20 * 60 -1)'; 

            %--- Combine the 20 values of each bit ------------------------ 
            bits = reshape(bits, 20, (size(bits, 1) / 20)); 
            bits = sum(bits); 

            % Now threshold and make it -1 and +1  
            bits(bits > 0)  = 1; 
            bits(bits <= 0) = -1; 

            %--- Check the parity of the TLM and HOW words ---------------- 
            if (navPartyChk(bits(1:32)) ~= 0) && ... 
               (navPartyChk(bits(31:62)) ~= 0) 
                % Parity was OK. Record the preamble start position. Skip 
                % the rest of preamble pattern checking for this channel 
                % and process next channel.  

                firstSubFrame(channelNr) = index(i); 
                break;     
            end % if parity is OK ... 

        end % if (~isempty(find(index2 == 6000))) 
    end % for i = 1:size(index) 
 
    % Exclude channel from the active channel list if no valid preamble was 
    % detected 
    if firstSubFrame(channelNr) == 0 
         
        % Exclude channel from further processing. It does not contain any 
        % valid preamble and therefore nothing more can be done for it. 
        activeChnList = setdiff(activeChnList, channelNr); 
         
        disp(['Could not find valid preambles in channel ', ... 
                                                  num2str(channelNr),'!']); 
    end 
     
end % for channelNr = activeChnList 


%% Decode Ephemerides
% 
fprintf('Decoding Ephemeris (Akos) \n')
 
for channelNr = activeChnList'
 
  %=== Convert tracking output to navigation bits ======================= 

  %--- Copy 5 sub-frames long record from tracking output --------------- 
  navBitsSamples = trackRes(channelNr).IP(firstSubFrame(channelNr) - 20 : ... 
                             firstSubFrame(channelNr) + (1500 * 20) -1)'; 

  %--- Group every 20 vales of bits into columns ------------------------ 
  navBitsSamples = reshape(navBitsSamples, ... 
                           20, (size(navBitsSamples, 1) / 20)); 

  %--- Sum all samples in the bits to get the best estimate ------------- 
  navBits = sum(navBitsSamples); 

  %--- Now threshold and make 1 and 0 ----------------------------------- 
  % The expression (navBits > 0) returns an array with elements set to 1 
  % if the condition is met and set to 0 if it is not met. 
  navBits = (navBits > 0); 

  %--- Convert from decimal to binary ----------------------------------- 
  % The function ephemeris expects input in binary form. In Matlab it is 
  % a string array containing only "0" and "1" characters. 
  navBitsBin = dec2bin(navBits); 

  %=== Decode ephemerides and TOW of the first sub-frame ================ 
  [eph(acq.svs(channelNr)), TOW] = ... 
                          ephemeris(navBitsBin(2:1501)', navBitsBin(1)); 

  %--- Exclude satellite if it does not have the necessary nav data ----- 
  if (isempty(eph(acq.svs(channelNr)).IODC) || ... 
    isempty(eph(acq.svs(channelNr)).IODE_sf2) || ... 
    isempty(eph(acq.svs(channelNr)).IODE_sf3)) 

    %--- Exclude channel from the list (from further processing) ------ 
    activeChnList = setdiff(activeChnList, channelNr); 
  end     
end 


%% Save which GPS times correspond 



%% End Matters
fprintf('Saving Results to `part2a_ephem.mat`\n')
save part2a_ephem

