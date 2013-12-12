%Acquisition
%Grant Apperson
 
clear all; close all; clc;
 
try
    if matlabpool('size') == 0 % checking to see if my pool is already open
        matlabpool open
    end
catch
end
 
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                               Figure Setup
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fig_num=8;
fig_rows=2;
fig_pos=get_fig_pos(fig_num,fig_rows);
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                               Constants
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
freqL1 = 154*10.23e6;
sample_frequency = 16.3676e6;
ts=1/sample_frequency;
intermediate_frequency = 4.1304e6;
integration_period = 1e-3;
chiprate=1023/1e-3; 
 
total_time=15; %sec
 
samp_per_chip=sample_frequency/chiprate;
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                               Load Data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% open data file
filename = '/Users/grantapperson/Desktop/Kingston/School/Classes/Grad School/Fall 2013/GPS/Class_GPS_data/run_north.sim';
 
 
% open file
fid = fopen(sprintf('%s',filename));
% number of bytes
bytes_to_read = ceil(sample_frequency*integration_period);
 
signal= fread(fid,bytes_to_read,'int8')';
 
fclose(fid);
 
num_samp=bytes_to_read;
dtau=.5;
chip_shift=round(dtau*samp_per_chip);
tau=(0:chip_shift:num_samp-1*chip_shift);
time=0:ts:ts*(num_samp-1);
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                   Generate CA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
CA_code=cacode([1:32],ceil(samp_per_chip));
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                           Determining Doppler freq
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%acquisition will be done from rest
dopp_max=5000; %5 kHz
dopp_min=-5000; %-5 kHz
num_dopp_bins=25;
diff_dopp=(dopp_max-dopp_min)/(num_dopp_bins-1);
dopp_freq=dopp_min:diff_dopp:dopp_max;
 
 
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                        Acquiring SVs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
out=zeros(length(tau),num_dopp_bins);
out_ratio = zeros(1,32);
for SV=1:32

for n=1:num_dopp_bins
    sine=imag(exp(1j*2*pi*(dopp_freq(n)+intermediate_frequency)*time));
    cosine=real(exp(1j*2*pi*(dopp_freq(n)+intermediate_frequency)*time));
     parfor k=1:length(tau)
         CA_shift=shift(CA_code(SV,:),tau(k));
         I=signal.*CA_shift.*sine;
         Q=signal.*CA_shift.*cosine;
         out(k,n)=sum(I)^2+sum(Q)^2;
     end
end
dopp_corr=max(out,[],1);
chip_corr=max(out,[],2);
[max_dopp(SV),best_dopp_bin]=max(dopp_corr);
[max_chip(SV),best_chip_shift]=max(chip_corr);
 
ave_out(SV)=mean(mean(out));
 
dopp_freq_aq=(best_dopp_bin-1)*diff_dopp+dopp_min;
chip_shift_aq=(best_chip_shift-1)*dtau;
 
Low_Res_Acquisition(1,SV)=dopp_freq_aq;
Low_Res_Acquisition(2,SV)=chip_shift_aq;
clear out
 
end
aq_out.dopp_freq=Low_Res_Acquisition(1,:);
aq_out.chip_shift=Low_Res_Acquisition(2,:);
aq_out.ave_out=ave_out;
aq_out.max_chip=max_chip;
 
save('run_north_low_res_acq.mat','aq_out')
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                               Find Svs present
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
 
SV_pres_thresh=12;
for SV=1:32
 
out_ratio(SV)=max_chip(SV)/ave_out(SV);
 
end
 
[out_ratio_,svs] = sort(out_ratio,'descend');
 
 
 
present_svs =svs(find(out_ratio_>SV_pres_thresh));
num_svs_pres=length(present_svs);
 
reorder_Low_Res_Acquisition=Low_Res_Acquisition(:,present_svs);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                           New Doppler freq
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
num_dopp_bins_new=40;
dopp_max_new=reorder_Low_Res_Acquisition(1,:)+ diff_dopp; %5 kHz
dopp_min_new=reorder_Low_Res_Acquisition(1,:)- diff_dopp; %-5 kHz
 
 
 
Gold_Code=cacode(present_svs,ceil(samp_per_chip));
num_samp=length(Gold_Code);
time=0:ts:ts*(num_samp-1);
 
for i=1:num_svs_pres
    
    diff_dopp_new(i)=(dopp_max_new(i)-dopp_min_new(i))/(num_dopp_bins_new-1);
    dopp_freq_new(i,:)=dopp_min_new(i):diff_dopp_new(i):dopp_max_new(i);
    CA_shift_corr(i,:)=shift(Gold_Code(i,:),reorder_Low_Res_Acquisition(2,i)*ceil(samp_per_chip));
    
    for n=1:num_dopp_bins_new
    sine=imag(exp(1j*2*pi*(dopp_freq_new(i,n)+intermediate_frequency)*time));
    cosine=real(exp(1j*2*pi*(dopp_freq_new(i,n)+intermediate_frequency)*time));
    I=signal.*CA_shift_corr(i,:).*sine;
    Q=signal.*CA_shift_corr(i,:).*cosine;
    out_fine(n)=sum(I)^2+sum(Q)^2;
    end
    
    dopp_corr_fine=max(out_fine,[],1);
    [max_dopp,best_dopp_bin]=max(dopp_corr_fine);
    dopp_freq_aq_fine(i)=(best_dopp_bin-1)*diff_dopp_new(i)+dopp_min_new(i);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                           save new Acquisition
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

run_north_high_res_acq.fdopp=dopp_freq_aq_fine;
run_north_high_res_acq.nsv=num_svs_pres;
run_north_high_res_acq.svs=present_svs;
run_north_high_res_acq.tau_chips=reorder_Low_Res_Acquisition(2,:);
run_north_high_res_acq.tau_samples=reorder_Low_Res_Acquisition(2,:)*ceil(samp_per_chip);

save('run_north_high_res_acq')
 


 
 
 





