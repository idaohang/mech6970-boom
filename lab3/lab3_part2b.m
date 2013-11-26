%Grant Apperson
%lab 3
%part 2b
%DGPS
%Relative Positioning

clear all; clc; close all;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                              Constants
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fig_num=8;
fig_rows=2;
fig_pos=get_fig_pos(fig_num,fig_rows);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                              Constants
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c=299792458; %m/s
L1_freq=1575.42*10^6; %Hz
L1_wv_len=c/L1_freq;  %m
rcvr_name={'gNovatel0','gNovatel1'};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                              Load Data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load(['..' filesep 'data' filesep 'GPS.mat'])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                    Removing SVs with partial data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
num_of_svs=length(GPS.(rcvr_name{1}).SV_Number);
avail_sv=GPS.(rcvr_name{1}).SV_Number([1:1,3:num_of_svs]);
    num_of_svs=length(avail_sv);
    for k=1:num_of_svs
        sv_num_str=num2str(avail_sv(k));

        sv_name=['SV' sv_num_str];
        if GPS.(rcvr_name{1}).PsrL1.(sv_name)(2)==0 || GPS.(rcvr_name{2}).PsrL1.(sv_name)(1)==0
            unavail_sv=k;
        end
        
    end
    avail_sv=avail_sv([1:unavail_sv-1,unavail_sv+1:num_of_svs]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                    Finding Position Solution
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
for count=1:length(GPS.(rcvr_name{2}).Time)
    count_rcvr1=count+1;

    num_of_svs_now=length(avail_sv);
    delta_param=ones(3,1);
    %---Intial Estimates (Center of Earth)---%
    user_pos_est=[0,0,0];
    b_est=0;

    
    pseudorange_est=zeros(num_of_svs_now,1);

    sv_pos_mtrx=zeros(num_of_svs_now,3);
 
    
    for k=1:num_of_svs_now
        sv_num_str=num2str(avail_sv(k));
        sv_name=['SV' sv_num_str];
        sv_pos_mtrx(k,:)=[GPS.(rcvr_name{1}).X_SV.(sv_name)(count_rcvr1),GPS.(rcvr_name{1}).Y_SV.(sv_name)(count_rcvr1),GPS.(rcvr_name{1}).Z_SV.(sv_name)(count_rcvr1)];
        
        
        pseudorange_est(k,1)=norm(sv_pos_mtrx(k,:)-user_pos_est)+b_est;
        %---Calculating Intial Pseudorange Errors---%
        delta_pseudorange(k,1)=(GPS.(rcvr_name{1}).PsrL1.(sv_name)(count_rcvr1)-pseudorange_est(k,1));
    end
    while norm(delta_param(1:3))>.001
        %---Calculating Geometry Matrix---%
        G=calc_geo_matrix(sv_pos_mtrx,user_pos_est);

        %---Calculating change in Estimates ---%
        delta_param=inv(G'*G)*G'*delta_pseudorange;
        
        %---Calculating New User Position and Clock Bias Estimates---%
        user_pos_est=user_pos_est+delta_param(1:3)';
        b_est=b_est+delta_param(4);

        %---Calculating New Pseudorange Estimates for Each SV---%
        for k=1:num_of_svs_now
            
            sv_num_str=num2str(avail_sv(k));
            sv_name=['SV' sv_num_str];
            
            pseudorange_est(k,1)=norm(sv_pos_mtrx(k,:)-user_pos_est)+b_est;
            rel_pseudorange(k,1)=GPS.(rcvr_name{1}).PsrL1.(sv_name)(count_rcvr1)-GPS.(rcvr_name{2}).PsrL1.(sv_name)(count);
            %---Calculating New Pseudorange Errors---%;
            delta_pseudorange(k,1)=(GPS.(rcvr_name{1}).PsrL1.(sv_name)(count_rcvr1)-pseudorange_est(k,1));
        end
    end
        del_rho(count,:)=pinv(G)*rel_pseudorange;
        range(count)=norm(del_rho(count,[1:3]));
        range_error(count)=1.905-range(count);

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                    Finding Position Solution
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for j=1:2
for count=1:length(GPS.(rcvr_name{j}).Time)
    for k=1:num_of_svs
        sv_num_str=num2str(GPS.(rcvr_name{j}).SV_Number(k));
        sv_name=['SV' sv_num_str];
        if GPS.(rcvr_name{j}).PsrL1.(sv_name)(count)==0;
            unavail_sv=k;
        end
    end
    avail_sv=GPS.(rcvr_name{j}).SV_Number([1:unavail_sv-1,unavail_sv+1:num_of_svs]);
    num_of_svs_now=length(avail_sv);
    delta_param=ones(3,1);
    %---Intial Estimates (Center of Earth)---%
    user_pos_est=[0,0,0];
    user_vel_est=[0,0,0];
    b_est=0;
    b_dot_est=0;
    
    pseudorange_est=zeros(num_of_svs_now,1);
    pseudorangerate=zeros(num_of_svs_now,1);
    sv_pos_mtrx=zeros(num_of_svs_now,3);
    sv_vel_mtrx=zeros(num_of_svs_now,3);
    
    for k=1:num_of_svs_now
        sv_num_str=num2str(avail_sv(k));
        sv_name=['SV' sv_num_str];
        sv_pos_mtrx(k,:)=[GPS.(rcvr_name{j}).X_SV.(sv_name)(count),GPS.(rcvr_name{j}).Y_SV.(sv_name)(count),GPS.(rcvr_name{j}).Z_SV.(sv_name)(count)];
        sv_vel_mtrx(k,:)=[GPS.(rcvr_name{j}).X_Vel_SV.(sv_name)(count),GPS.(rcvr_name{j}).Y_Vel_SV.(sv_name)(count),GPS.(rcvr_name{j}).Z_Vel_SV.(sv_name)(count)];
        
        pseudorange_est(k,1)=norm(sv_pos_mtrx(k,:)-user_pos_est)+b_est;
        pseudorangerate(k,1)=GPS.(rcvr_name{j}).DopL1.(sv_name)(count)*L1_wv_len;
        psr(k,1)=GPS.(rcvr_name{j}).PsrL1.(sv_name)(count);
        %---Calculating Intial Pseudorange Errors---%
        delta_pseudorange(k,1)=(GPS.(rcvr_name{j}).PsrL1.(sv_name)(count)-pseudorange_est(k,1));
    end
    while norm(delta_param(1:3))>.001
        %---Calculating Geometry Matrix---%
        G=calc_geo_matrix(sv_pos_mtrx,user_pos_est);
        unit_vect=G(:,1:3)';
        
        for k=1:num_of_svs_now
            pseudorangerate_tilde(k,1)=pseudorangerate(k,1)-sv_vel_mtrx(k,:)*unit_vect(:,k);
        end
        %---Calculating change in Estimates ---%
        delta_param=inv(G'*G)*G'*delta_pseudorange;
        vel_param=inv(G'*G)*G'*pseudorangerate_tilde;
        
        %---Calculating New User Position and Clock Bias Estimates---%
        user_pos_est=user_pos_est+delta_param(1:3)';
        b_est=b_est+delta_param(4);
        
        user_vel_est=vel_param(1:3)';
        b_dot_est=vel_param(4);
        
        %---Calculating New Pseudorange Estimates for Each SV---%
        for k=1:num_of_svs_now
            
            sv_num_str=num2str(avail_sv(k));
            sv_name=['SV' sv_num_str];
            
            pseudorange_est(k,1)=norm(sv_pos_mtrx(k,:)-user_pos_est)+b_est;
        
            %---Calculating New Pseudorange Errors---%;
            delta_pseudorange(k,1)=(GPS.(rcvr_name{j}).PsrL1.(sv_name)(count)-pseudorange_est(k,1));
        end
    end
    GPS.(rcvr_name{j}).X_Rcvr(count)=user_pos_est(1);
    GPS.(rcvr_name{j}).Y_Rcvr(count)=user_pos_est(2);
    GPS.(rcvr_name{j}).Z_Rcvr(count)=user_pos_est(3);
    GPS.(rcvr_name{j}).X_Vel_Rcvr(count)=user_vel_est(1);
    GPS.(rcvr_name{j}).Y_Vel_Rcvr(count)=user_vel_est(2);
    GPS.(rcvr_name{j}).Z_Vel_Rcvr(count)=user_vel_est(3);
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                    Convert Position Solution to LLA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [GPS.(rcvr_name{j}).Lat_Rcvr(count),GPS.(rcvr_name{j}).Lon_Rcvr(count),GPS.(rcvr_name{j}).Alt_Rcvr(count)]=Wgsxyz2lla(user_pos_est');
    [GPS.(rcvr_name{j}).Lat_Meas(count),GPS.(rcvr_name{j}).Lon_Meas(count),GPS.(rcvr_name{j}).Alt_Meas(count)]=Wgsxyz2lla([GPS.(rcvr_name{j}).X_Meas(count),GPS.(rcvr_name{j}).Y_Meas(count),GPS.(rcvr_name{j}).Z_Meas(count)]);
    ENU_Rcvr=Wgsxyz2enu(user_pos_est',GPS.(rcvr_name{j}).Lat_Meas(count),GPS.(rcvr_name{j}).Lon_Meas(count),GPS.(rcvr_name{j}).Alt_Meas(count));
    GPS.(rcvr_name{j}).East_Rcvr(count)=ENU_Rcvr(1);
    GPS.(rcvr_name{j}).North_Rcvr(count)=ENU_Rcvr(2);
    GPS.(rcvr_name{j}).Up_Rcvr(count)=ENU_Rcvr(3);
end
end

for count=1:length(GPS.(rcvr_name{2}).Time)
    count_rcvr1=count+1;
   
    pos_diff(count)=norm([GPS.(rcvr_name{1}).X_Rcvr(count_rcvr1)-GPS.(rcvr_name{2}).X_Rcvr(count),GPS.(rcvr_name{1}).Y_Rcvr(count_rcvr1)-GPS.(rcvr_name{2}).Y_Rcvr(count),GPS.(rcvr_name{1}).Z_Rcvr(count_rcvr1)-GPS.(rcvr_name{2}).Z_Rcvr(count)]);
    pos_diff_error(count)=1.905-pos_diff(count);
end

pos_diff_error_mean=mean(pos_diff_error);
pos_diff_error_std=std(pos_diff_error);

range_error_mean=mean(range_error);
range_error_std=std(range_error);
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                               Plotting
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fig(1)=figure;set(fig(1),'Position',fig_pos{1})
hold on
title('Range Between Antennas')
plot(GPS.(rcvr_name{2}).Time,pos_diff)
plot(GPS.(rcvr_name{2}).Time,range,'--r')
xlabel('Time (sec)')
ylabel('Range (m)')
legend('Receiver Difference','DGPS Solution','Location','Best')
hold off

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                               Printing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('Part2 a) \n')
fprintf('Mean Error of Difference in Receiver Position: %f (m)\n',pos_diff_error_mean)
fprintf('Standard Deviation of Error in Difference Between Receiver Positions: %f (m)\n',pos_diff_error_std)
fprintf('-----------------------------\n')
fprintf('Part2 b) \n')
fprintf('Mean Error of DGPS Position Solution: %f \n',range_error_mean)
fprintf('Standard Deviation of Error in DGPS Position Solution: %f \n',range_error_std)

