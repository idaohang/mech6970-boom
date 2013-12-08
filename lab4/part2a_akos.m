%% MECH 6970 Lab4, Part 2, (a) - Primary, using DM Akos' code for comparison
% 
% Loads the results of the acquisition refinement and does tracking
% 
% @author Robert Cofield
% 
tic
clear all; clc
acq = load('part2a_narrow_ack.mat');

filename = ['..' filesep 'data' filesep 'GPS_Data_NordNav1e.sim'];
fileid = fopen(filename);

addpath(['..' filesep 'akos']);


%% Constants

svidx_ = 1;

fL1 = 154*10.23e6; % L1 frequency,  1.5754e+09 Hz
fs = 16.3676e6; % sampling frequency
fIF = 4.1304e6; % intermediate frequency
fD = 50; % Data message chipping frequency, Hz

TL1 = 1/fL1; % L1 carrier period
Tca_chip = 1.023e-6; % L1 C/A chip period
Tca = 1.0e-3; % L1 C/A sequence period
Td_chip = 1/fD;
Ts = 1/fs; % sampling period
Tid = Tca; % integrate & dump period (sec) - use the C/A code period
Tnav_chip = 20e-3; % period of a single nav msg chip
Tnav_frame = Tnav_chip*1500; % period of the entire nav msg frame
Tnav_subframe = Tnav_frame/5; % period of a nav msg subframe

stop_time = 1; % how long to go in sec

preamble = bin2dec('10001011');  

% % DLL filter 2nd order
% K_dll = 0.01;
% a_dll = -0.75;

% PLL filter 2nd Order
pll = struct(...
  'z',0.9,...
  'f',5);

% % Look at PLL filter
% num = [4*pi*pll.z*pll.f , (2*pi*pll.f)^2];
% den = [1 , 4*pi*pll.z*pll.f ,  (2*pi*pll.f)^2];
% rlocus(tf(num,den,Tid))
% pause

