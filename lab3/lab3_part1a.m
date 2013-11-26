%Grant Apperson
%lab 3
%part 1c
%Velocity Estimation

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
%data=load('/Volumes/KINGSTON/School/Classes/Grad School/Fall 2013/GPS/Class_GPS_data/Novatel_parsed_data_2.mat');
data=load(['..' filesep 'data' filesep 'Novatel_Data_ephemfixed.mat']);

for j=1:length(rcvr_name);
    %+++++++++++++++++Creating Strings to Parse Structs++++++++++++++++++++++%
    
    %-----Raw Data-----%
    psrL1_struct_string=['zPsrL1_' rcvr_name{j}];
    dopL1_struct_string=['zDopL1_' rcvr_name{j}];
    adrL1_struct_string=['zAdrL1_' rcvr_name{j}];
    cnoL1_struct_string=['zCnoL1_' rcvr_name{j}];
    
    %-----Position Solution of Receiver-----%
    zX_struct_string=['zX_' rcvr_name{j}];
    zY_struct_string=['zY_' rcvr_name{j}];
    zZ_struct_string=['zZ_' rcvr_name{j}];
    
    %-----Velocity Solution of Receiver-----%
    zVelX_struct_string=['zVelX_' rcvr_name{j}];
    zVelY_struct_string=['zVelY_' rcvr_name{j}];
    zVelZ_struct_string=['zVelZ_' rcvr_name{j}];
    
    %+++++++++++++++++++++++++Parse Structs+++++++++++++++++++++++++++++++%
    
    %-----Raw Data-----%
    psrL1_=data.(rcvr_name{j}).(psrL1_struct_string);
    dopL1_=data.(rcvr_name{j}).(dopL1_struct_string);
    adrL1_=data.(rcvr_name{j}).(adrL1_struct_string);
    cnoL1_=data.(rcvr_name{j}).(cnoL1_struct_string);
    
    %-----Position Solution of Receiver-----%
    X_meas=data.(rcvr_name{j}).(zX_struct_string);
    Y_meas=data.(rcvr_name{j}).(zY_struct_string);
    Z_meas=data.(rcvr_name{j}).(zZ_struct_string);
    
    %-----Velocity Solution of Receiver-----%
    X_vel_meas=data.(rcvr_name{j}).(zVelX_struct_string);
    Y_vel_meas=data.(rcvr_name{j}).(zVelY_struct_string);
    Z_vel_meas=data.(rcvr_name{j}).(zVelZ_struct_string);
    
    
    %++++++++++Find Samples Without all required Measurements+++++++++++++%
    num_epoch=length(X_meas);
    act_data = [];
    
    %-----Find Rows That Have all Data-----%
    for count=1:num_epoch
        if isempty(psrL1_{count}) || isempty(dopL1_ {count}) || isnan(X_meas (count))
            continue
        end
        time_epoch(count) = psrL1_{count}(34); %time measurement (time of week (ms) ) for each pseudorange measurement
        act_data(end+1) = count;
    end
    
    %-----Update Raw Data-----%
    time_epoch=time_epoch(act_data);
    psrL1_=psrL1_(act_data);
    dopL1_=dopL1_(act_data);
    adrL1_=adrL1_(act_data);
    cnoL1_=cnoL1_(act_data);
    
    %-----Update Position Solution-----%
    X_meas=X_meas(act_data);
    Y_meas=Y_meas(act_data);
    Z_meas=Z_meas(act_data);
    
    %-----Update Velocity Solution-----%
    X_vel_meas=X_vel_meas(act_data);
    Y_vel_meas=Y_vel_meas(act_data);
    Z_vel_meas=Z_vel_meas(act_data);
    
    num_epoch=length(time_epoch);
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                           Find SVs in View
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %+++++++++++++++++++++++++++++Try all SVs ++++++++++++++++++++++++++++%
    count=1;
    for k=1:32
        sv_num_str=num2str(k);
        sv_name=['SV' sv_num_str];
        data_name=['zEphem' sv_num_str '_' rcvr_name{j}];
        try
            
            Ephem.(sv_name)= data.(rcvr_name{j}).(data_name){end};
            sv_num(count)=k;
            count=count+1;
        catch
            Ephem.(sv_name)=zeros(30,1);
        end
        
    end
    
    num_of_svs=length(sv_num);
    
    %+++++++++++++++++++Transfer RAW Data to Matrix ++++++++++++++++++++++%
    
    %-----Preallocating-----%
    psrL1=zeros(num_of_svs,num_epoch);
    dopL1=zeros(num_of_svs,num_epoch);
    adrL1=zeros(num_of_svs,num_epoch);
    cnoL1=zeros(num_of_svs,num_epoch);
    
    %-----Rearranging-----%
    for count= 1:num_epoch
        psrL1(:,count)=psrL1_{count}(sv_num+1); 
        dopL1(:,count)=dopL1_{count}(sv_num+1);
        adrL1(:,count)=adrL1_{count}(sv_num+1);
        cnoL1(:,count)=cnoL1_{count}(sv_num+1);
           
    end
    
    %+++++++++++++++++Only Use Raw Data With Ephem and Psr Data ++++++++++++++++++%
    have_dat = find(any(psrL1,2));
    sv_num=sv_num(have_dat);
    psrL1=psrL1(have_dat,:);
    dopL1=dopL1(have_dat,:);
    adrL1=adrL1(have_dat,:);
   
    
    num_of_svs=length(sv_num);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                              Create GPS Struct
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %---Preallocating SV Specific Data---%
    GPS.(rcvr_name{j}).Ephem=[];
    GPS.(rcvr_name{j}).PsrL1=[];
    GPS.(rcvr_name{j}).DopL1=[];
    GPS.(rcvr_name{j}).AdrL1=[];
    GPS.(rcvr_name{j}).CnoL1=[];
    GPS.(rcvr_name{j}).X_SV=[];
    GPS.(rcvr_name{j}).Y_SV=[];
    GPS.(rcvr_name{j}).Z_SV=[];
    GPS.(rcvr_name{j}).X_Vel_SV=[];
    GPS.(rcvr_name{j}).Y_Vel_SV=[];
    GPS.(rcvr_name{j}).Z_Vel_SV=[];
    
    %---Preallocating Data to be Calculated---%
    GPS.(rcvr_name{j}).X_Rcvr=[];
    GPS.(rcvr_name{j}).Y_Rcvr=[];
    GPS.(rcvr_name{j}).Z_Rcvr=[];
    
    GPS.(rcvr_name{j}).Lat_Rcvr=[];
    GPS.(rcvr_name{j}).Lon_Rcvr=[];
    GPS.(rcvr_name{j}).Alt_Rcvr=[];
    
    GPS.(rcvr_name{j}).X_Vel_Rcvr=[];
    GPS.(rcvr_name{j}).Y_Vel_Rcvr=[];
    GPS.(rcvr_name{j}).Y_Vel_Rcvr=[];
    
    
    %++++++++++++++++++ Storing Initial Data in Struct +++++++++++++++++++%
    GPS.(rcvr_name{j}).SV_Number=sv_num;    
    GPS.(rcvr_name{j}).Time=time_epoch/1000;
    GPS.(rcvr_name{j}).X_Meas=X_meas;
    GPS.(rcvr_name{j}).Y_Meas=Y_meas;
    GPS.(rcvr_name{j}).Z_Meas=Z_meas;
    GPS.(rcvr_name{j}).X_Vel_Meas=X_vel_meas;
    GPS.(rcvr_name{j}).Y_Vel_Meas=Y_vel_meas;
    GPS.(rcvr_name{j}).Z_Vel_Meas=Z_vel_meas;
    
    
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                       Caclulate SVs Position and Velocity
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   for k=1:num_of_svs
        sv_num_str=num2str(sv_num(k));
        sv_name=['SV' sv_num_str];
        
        for count=1:num_epoch
            transitTime=psrL1(k,count)/c;
            transmitTime=time_epoch(count)/1000-transitTime;
            
            
            if psrL1(k,count)==0
                GPS.(rcvr_name{j}).X_SV.(sv_name)(count)=NaN;
                GPS.(rcvr_name{j}).Y_SV.(sv_name)(count)=NaN;
                GPS.(rcvr_name{j}).Z_SV.(sv_name)(count)=NaN;
                GPS.(rcvr_name{j}).X_Vel_SV.(sv_name)(count)=NaN;
                GPS.(rcvr_name{j}).Y_Vel_SV.(sv_name)(count)=NaN;
                GPS.(rcvr_name{j}).Z_Vel_SV.(sv_name)(count)=NaN;
                continue
            end
            
            [SV_Pos,SV_Vel,SV_Clk_corr]=calc_sv_pos_and_vel(Ephem.(sv_name),transmitTime,transitTime);
            GPS.(rcvr_name{j}).X_SV.(sv_name)(count)=SV_Pos(1);
            GPS.(rcvr_name{j}).Y_SV.(sv_name)(count)=SV_Pos(2);
            GPS.(rcvr_name{j}).Z_SV.(sv_name)(count)=SV_Pos(3);
            GPS.(rcvr_name{j}).X_Vel_SV.(sv_name)(count)=SV_Vel(1);
            GPS.(rcvr_name{j}).Y_Vel_SV.(sv_name)(count)=SV_Vel(2);
            GPS.(rcvr_name{j}).Z_Vel_SV.(sv_name)(count)=SV_Vel(3);
            
            GPS.(rcvr_name{j}).PsrL1.(sv_name)(count)=psrL1(k,count)+SV_Clk_corr*c;
        end
    
        GPS.(rcvr_name{j}).Ephem.(sv_name)=Ephem.(sv_name);
        GPS.(rcvr_name{j}).DopL1.(sv_name)=dopL1(k,:);
        GPS.(rcvr_name{j}).AdrL1.(sv_name)=adrL1(k,:);
        GPS.(rcvr_name{j}).CnoL1.(sv_name)=cnoL1(k,:);
        
    end

    clear psrL1_ dopL1_ adrL1_ cnoL1_ time_epoch act_data num_epoch
end


j=1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                   Plotting Constellation of SVs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fig(1)=figure;set(fig(1),'Position',fig_pos{1})
hold on
for k=1:num_of_svs
    sv_num_str=num2str(GPS.(rcvr_name{j}).SV_Number(k));
    sv_name=['SV' sv_num_str];
    scatter3(GPS.(rcvr_name{j}).X_SV.(sv_name),GPS.(rcvr_name{j}).Y_SV.(sv_name),GPS.(rcvr_name{j}).Z_SV.(sv_name),'r')
    
end
scatter3(0,0,0)

r = 6371e3; %earth's rad
[x_earth,y_earth,z_earth] = sphere(50);
x_earth = x_earth*r;
y_earth = y_earth*r;
z_earth = z_earth*r;

lightGrey = 0.8*[1 1 1];
surface(x_earth,y_earth,z_earth,'FaceColor', 'none','EdgeColor',lightGrey)


hold off
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                       clearing Unused Variables
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clearvars -except c L1_freq L1_wv_len rcvr_name GPS fig_pos
j=1; %only using gNovatel0
num_of_svs=length(GPS.(rcvr_name{j}).SV_Number);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                    Finding Position Solution
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                               Plotting
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%++++++++++++++++++++++++++++++ECEF Position++++++++++++++++++++++++++++++%
fig(2)=figure;set(fig(2),'Position',fig_pos{2})
%--------------X------------%
subplot(3,1,1)
hold on
plot(GPS.(rcvr_name{j}).Time,GPS.(rcvr_name{j}).X_Rcvr,'-b')
plot(GPS.(rcvr_name{j}).Time,GPS.(rcvr_name{j}).X_Meas,'--r')
xlabel('Time')
ylabel('ECEF X (m)')
legend('Calculated Value','Measured from Novatel')
hold off
%--------------Y------------%
subplot(3,1,2)
hold on
plot(GPS.(rcvr_name{j}).Time,GPS.(rcvr_name{j}).Y_Rcvr,'-b')
plot(GPS.(rcvr_name{j}).Time,GPS.(rcvr_name{j}).Y_Meas,'--r')
xlabel('Time')
ylabel('ECEF Y (m)')
legend('Calculated Value','Measured from Novatel')
hold off
%--------------Z------------%
subplot(3,1,3)
hold on
plot(GPS.(rcvr_name{j}).Time,GPS.(rcvr_name{j}).Z_Rcvr,'-b')
plot(GPS.(rcvr_name{j}).Time,GPS.(rcvr_name{j}).Z_Meas,'--r')
xlabel('Time')
ylabel('ECEF Z (m)')
legend('Calculated Value','Measured from Novatel')
hold off

%+++++++++++++++++++++++++ECEF Position Errors++++++++++++++++++++++++++++%

fig(3)=figure;set(fig(3),'Position',fig_pos{3})
%--------------X------------%
subplot(3,1,1)
hold on
plot(GPS.(rcvr_name{j}).Time,GPS.(rcvr_name{j}).X_Meas'-GPS.(rcvr_name{j}).X_Rcvr,'-b')
xlabel('Time')
ylabel('ECEF X Error (m)')
hold off
%--------------Y------------%
subplot(3,1,2)
hold on
plot(GPS.(rcvr_name{j}).Time,GPS.(rcvr_name{j}).Y_Meas'-GPS.(rcvr_name{j}).Y_Rcvr,'-b')
xlabel('Time')
ylabel('ECEF Y Error (m)')
hold off
%--------------Z------------%
subplot(3,1,3)
hold on
plot(GPS.(rcvr_name{j}).Time,GPS.(rcvr_name{j}).Z_Meas'-GPS.(rcvr_name{j}).Z_Rcvr,'-b')
xlabel('Time')
ylabel('ECEF Z Error(m)')
hold off


%+++++++++++++++++++++++++ENU Position Errors++++++++++++++++++++++++++++%

fig(4)=figure;set(fig(4),'Position',fig_pos{4})
%--------------East------------%
subplot(3,1,1)
hold on
plot(GPS.(rcvr_name{j}).Time,GPS.(rcvr_name{j}).East_Rcvr,'-b')
xlabel('Time')
ylabel('East Error (m)')
hold off
%--------------North------------%
subplot(3,1,2)
hold on
plot(GPS.(rcvr_name{j}).Time,GPS.(rcvr_name{j}).North_Rcvr,'-b')
xlabel('Time')
ylabel('North Error (m)')
hold off
%--------------Up------------%
subplot(3,1,3)
hold on
plot(GPS.(rcvr_name{j}).Time,GPS.(rcvr_name{j}).Up_Rcvr,'-b')
xlabel('Time')
ylabel('Up Error(m)')
hold off


%+++++++++++++++++++++++++Velocity Errors++++++++++++++++++++++++++++%

fig(5)=figure;set(fig(5),'Position',fig_pos{5})
%--------------X------------%
subplot(3,1,1)
hold on
plot(GPS.(rcvr_name{j}).Time,GPS.(rcvr_name{j}).X_Vel_Rcvr,'-b')
plot(GPS.(rcvr_name{j}).Time,GPS.(rcvr_name{j}).X_Vel_Meas,'--r')
xlabel('Time')
ylabel('ECEF X (m/s)')
legend('Calculated Value','Measured from Novatel')
hold off
%--------------Y------------%
subplot(3,1,2)
hold on
plot(GPS.(rcvr_name{j}).Time,GPS.(rcvr_name{j}).Y_Vel_Rcvr,'-b')
plot(GPS.(rcvr_name{j}).Time,GPS.(rcvr_name{j}).Y_Vel_Meas,'--r')
xlabel('Time')
ylabel('ECEF Y (m/s)')
legend('Calculated Value','Measured from Novatel')
hold off
%--------------Z------------%
subplot(3,1,3)
hold on
plot(GPS.(rcvr_name{j}).Time,GPS.(rcvr_name{j}).Z_Vel_Rcvr,'-b')
plot(GPS.(rcvr_name{j}).Time,GPS.(rcvr_name{j}).Z_Vel_Meas,'--r')
xlabel('Time')
ylabel('ECEF Z (m/s)')
legend('Calculated Value','Measured from Novatel')
hold off

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                               Printing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('-----------------------------\n')
fprintf('Part1 c) \n')
fprintf('Doppler measurements were used to find the velocity of the antenna\n because carrier measurements are much more accurate when \n calcualting a differenece.')  

