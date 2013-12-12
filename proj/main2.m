%% Numerically Integrate Crossbow

clear
load workspace_after_time_sync

% XBow_accY_filt = zeros(size(XBow_accY));
%h = fdesign.lowpass('Fp,Fst,Ap,Ast', 0.15, );
%d = design(h,'equiripple');
[b,a]=butter(4, .01, 'low');

for k = 1:4
  XBow_accY_filt(:,k) = filtfilt(b,a,XBow_accY(:,k));
end

plot(XBow_accY(:,1)); hold on
plot(XBow_accY_filt(:,1), 'g')


%% Our Tracking comparison
