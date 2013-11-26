%% MECH 6970 Lab4 Part 2, a) Preliminary 
% Narrow the Acquisition results of Part 1, a), for SV's which are determined to
% be present
% 
% This assumes you downloaded the data file into ../data/
% 
% NOTE !!! :
% The bit-level operations and array sizing for this 1ms analysis should not be
% exactly replicated for the 10ms part.
%   - in particular, the calculation of the upsample rate needs to be redone.
% 
% NOTE!!! : 
% This file extends Part 1a), so any improvements made to Part 1b) need to be
% implemented here as well
% 
% Loads the results of the serial-search-based acquisition from Part 1a),
% performed on all SV's
% 
% @author Robert Cofield
% 
tic 
clear all; close all; clc;
fprintf('Part II a) - Narrowing Acquisition Results\n')
% try matlabpool 3; catch e, disp(e); end

% Load the Acquisition data for all SV's as it was processed by Part I a)
load part1a_all_sv


%% Constants & Settingsfilename = ['..' filesep 'data' filesep 'GPS_Data_NordNav1e.sim'];

corr_thold_sv_present = 17;
tau_step = 1; % in actual measured samples
fdopp_step = 25; % Hz
sv_to_plot = 17; % which of the present sv's to use

%% Find SV's which are present
% Should be 1,2,4,8,9,10,12,17,20,24,28,32 ?
% By looking at crosscorrelation plots manually:
%   4,17,28,2,20

% % Use this to see which SV's should be present, choose threshold
% fh = figure;
%   plot(y_ratio)
%   grid on
%   title('Ratio of max correlation value to mean correlation value for each SV')
%   xlabel('SV #')
  
% % Use this to examine crosscorrelations
% fh = figure;
%   surf(fdopp,tau,y{20})
%   xlabel('Doppler Frequency (Hz)'); ylabel('Sample Shifts'); zlabel('Correlation');
%   title('CrossCorrelation of Signal with Replica');

% reorder with strongest correlation first
[y_ratio_, present_svs] = sort(y_ratio,'descend');
present_svs = present_svs(find(y_ratio_>corr_thold_sv_present));
nsv = length(present_svs);
clear y_ratio_


%% Hone in on Doppler Frequency and Arrival Time

yp = cell(1,nsv); % y for present sv
tau_p = cell(1,nsv); % narrowed range of taus used for each present sv
fdopp_p = cell(1,nsv); % narrowed range of fdopps used for each present sv
tau_p_soln = zeros(1,nsv);
tau_p_soln_idx = zeros(1,nsv);
fdopp_p_soln = zeros(1,nsv);
fdopp_p_soln_idx = zeros(1,nsv);
for svpi = 1:length(present_svs)
  sv = present_svs(svpi);
  
  % get a narrower range of tau values - between
  tau_lo = tau(tau_soln_idx(sv)-2);
  tau_hi = tau(tau_soln_idx(sv)+2);
  tau_ = tau_lo:tau_step:tau_hi;
  
  % get a narrower range of fdopps  
  feff_lo = feff(fdopp_soln_idx(sv)-2);
  feff_hi = feff(fdopp_soln_idx(sv)+2);
  feff_ = feff_lo:fdopp_step:feff_hi;
  
  % new landscape to calculate
  yp_ = zeros(length(tau_),length(feff_));
  
  % search over new tau range
  for t_ = 1:length(tau_)
    prn_shifted = shift(prn(sv,:),tau_(t_));
    % search over new fdopp (feff) range
    for f_ = 1:length(feff_); 
      sin_ = imag(exp(1j*2*pi*feff_(f_)*T));
      cos_ = real(exp(1j*2*pi*feff_(f_)*T));
      I = signal1.*prn_shifted.*sin_;
      Q = signal1.*prn_shifted.*cos_;
      yp_(t_,f_) = sum(I)^2 + sum(Q)^2;
    end    
  end
  
  % save data for this SV
  yp{1,svpi} = yp_;
  tau_p{1,svpi} = tau_;
  fdopp_p{1,svpi} = feff_ -fIF;
  
  % find the answer
%   [max_, fdopp_max_idx] = max(y{1,sv},[],2);
%   [~, tau_max_idx] = max(max_);
%   fdopp_soln_idx(sv) = fdopp_max_idx(tau_max_idx);
%   fdopp_soln(sv) = fdopp(fdopp_soln_idx(sv));
%   tau_soln_idx(sv) = tau_max_idx;
%   tau_soln(sv) = tau(tau_soln_idx(sv));
  [max_, fdopp_max_idx] = max(yp_,[],2);
  [~,tau_max_idx] = max(max_);
  fdopp_p_soln_idx(svpi) = fdopp_max_idx(tau_max_idx);
  fdopp_p_soln(svpi) = fdopp_p{1,svpi}(fdopp_p_soln_idx(svpi));
  tau_p_soln_idx(svpi) = tau_max_idx;
  tau_p_soln(svpi) = tau_(tau_max_idx);

end

clear feff_ tau_ yp_ tau_lo tau_hi feff_lo feff_hi I Q sin_ cos_ t_ f_
clear fdopp fdopp_bound fdopp_max_idx fdopp_soln fdopp_soln_idx fdopp_step
clear max_ feff tau_soln tau_soln_idx tau_step svpi y tau tau_chip_size tau_idx
clear tau_max_idx y_ratio y_max y_mean sv e nfdopp prn_shifted

%% Examine Isolated Doppler Frequency & Arrival Time

sv_to_plot_idx = find(present_svs==sv_to_plot);
fh = figure;
  surf(fdopp_p{sv_to_plot_idx},tau_p{sv_to_plot_idx},yp{sv_to_plot_idx})
  xlabel('Doppler Frequency (Hz)'); ylabel('Sample Shifts'); zlabel('Correlation');
  title(['Narrowed Search, SV #' num2str(sv_to_plot) ' - CrossCorrelation of Signal with Replica']);
  hold on
  plot3(fdopp_p_soln(sv_to_plot_idx), ...
        tau_p_soln(sv_to_plot_idx), ...
        yp{sv_to_plot_idx}(tau_p_soln_idx(sv_to_plot_idx),fdopp_p_soln_idx(sv_to_plot_idx)), ...
        'm.','MarkerSize',20 ...
      );
saveas(fh,'part2a_xcorr_closeup.png');
saveas(fh,'part2a_xcorr_closeup.fig');

clear sv_to_plot_idx fh sv_to_plot


%% End matters
% try matlabpool close; catch e, disp(e); end
% % Rename
tau_samples = tau_p_soln;
tau_chips = tau_samples/upsample;
fdopp = fdopp_p_soln;
svs = present_svs;

clear yp signal1 signal2 tau_p fdopp_p tau_p_soln_idx fdopp_p_soln_idx dfdopp
clear corr_thold_sv_present upsample integration_period fs prn N PC1 T Tca Ts
clear fIF fL1 tau_p_soln fdopp_p_soln present_svs

save part2a_narrow_ack
toc

