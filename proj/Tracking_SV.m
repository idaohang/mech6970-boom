%Tracking

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

total_time=120; %sec

samp_per_chip=sample_frequency/chiprate;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                               Load Data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load('run_north_high_res_acq.mat')

% open data file
filename = ['..' filesep 'data' filesep 'run_north.sim'];


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

dopp_freq(1)=run_north_high_res_acq.fdopp(1);
SV=run_north_high_res_acq.svs(1);
chip_shift_aq=run_north_high_res_acq.tau_chips(1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                           Tracking SVs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



dopp_freq(1)=560;
theta_error(1)=0;
theta(1)=0;
dopp_correction(1)=0;


K=.1;
a=.25;
Gold_code=shift(cacode(SV,ceil(samp_per_chip)),(chip_shift_aq)*ceil(samp_per_chip));

fid = fopen(sprintf('%s',filename));

    for k=1:total_time/integration_period
        signal= fread(fid,bytes_to_read,'int8')';
        
        %---Generate Shifted Gold Code---%
            CA_code_P=Gold_code;
            CA_code_E=shift(CA_code_P,chip_shift);
            CA_code_L=shift(CA_code_P,-chip_shift);
            
        %---NCO---%
            sine=imag(exp(1j*(2*pi*(dopp_freq(k)+intermediate_frequency).*time+theta(k))));
            cosine=real(exp(1j*(2*pi*(dopp_freq(k)+intermediate_frequency).*time +theta(k))));
            
        %---Delay Lock---%
            IE(k)=sum(signal.*sine.*CA_code_E);
            IL(k)=sum(signal.*sine.*CA_code_L);
            IP(k)=sum(signal.*sine.*CA_code_P);
            
            QE(k)=sum(signal.*cosine.*CA_code_E);
            QL(k)=sum(signal.*cosine.*CA_code_L);
            QP(k)=sum(signal.*cosine.*CA_code_P);
            
            %plot([-.5,0,.5],[(IE(k))^2+(QE(k))^2,(IP(k))^2+(QP(k))^2,(IL(k))^2+(QL(k))^2])
            
            code_error(k)=(sqrt(IE(k)^2+QE(k)^2)-sqrt(IL(k)^2+QL(k)^2))/(sqrt(IE(k)^2+QE(k)^2)+sqrt(IL(k)^2+QL(k)^2));
            
        %---Phase Lock---%
            Inphase(k)=sum(signal.*sine.*CA_code_P);
            Quadrature(k)=sum(signal.*cosine.*CA_code_P);
             
            theta_error(k+1)=atan(Quadrature(k)/Inphase(k));
            
        %generate new dopp freq from thete error 
            %dopp_correction(k+1)=dopp_correction(k)+K1*(theta_error(k+1)-theta_error(k))+K2*theta_error(k+1);
            dopp_correction(k+1)=dopp_correction(k)+K*(theta_error(k+1)-a*theta_error(k));
            dopp_freq(k+1)=dopp_freq(1)+dopp_correction(k+1);
            
            theta(k+1)=rem(2*pi*(dopp_freq(k+1)+intermediate_frequency)*(time(end)+ts)+theta(k),2*pi);
            
            shift_(k)=round(samp_per_chip*code_error(k));
            Gold_code=shift(Gold_code,shift_(k)); 
    end





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                               Plotting
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


fig1=figure;set(fig1,'Position',fig_pos{1})
hold on
plot([-.5,0,.5],[(IE(k))^2+(QE(k))^2,(IP(k))^2+(QP(k))^2,(IL(k))^2+(QL(k))^2])
hold off


fig2=figure;set(fig2,'Position',fig_pos{2})
hold on
title('In-Phase')
plot([0:integration_period:total_time-1*integration_period],Inphase)
xlabel('Time (sec)')
hold off

fig3=figure;set(fig3,'Position',fig_pos{3})
hold on
title('Doppler Freq. ')
plot([0:integration_period:total_time],dopp_freq)
xlabel('Time (sec)')
hold off
