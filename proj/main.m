%% main.m Primary File for Final Project
% 
% 
clear all; close all; clc
tic

%% Settings / Constants

% run = 'south';
run = 'north';

% filter_order = 1;
filter_order = 2;
% filter_order = 3;

data_dir = ['..' filesep 'data' filesep 'final_proj_data' filesep];


%% load data

load(['run_' run '_alog_interp.mat']);
nordnav_fid = fopen([data_dir 'run_' run '.sim'], 'r');


%% Get SV Geometry Data

% use the Novatel and find it for each GPS time epoch



%% Do Acquisition
% This needs to output a struct that contains the same data as `acq` in the file
% `../lab4/part2a_narrow_ack.m`




%% Do tracking with one of the filters
% This needs to output a struct that contains the same data as `trackRes` in the
% file `../lab4/part2a_tracking.m`

switch filter_order
  case 1
    % do 1st order filtering
  case 2
    % do 2nd order filtering
  case 3
    % do 3rd order filtering
end


%% Get GPS times from the data

% % number of milliseconds of data (CA code periods)
% n_code_per = length(trackRes(1).IP);
% for k = 1:length(trackRes)
%   if length(trackRes(k).IP) ~= n_code_per
%     error('Nonequal IP lengths');
%   end
% end
% 
% [TLM_starts, TOW, ch_status] = getTOWatIPIdx(trackRes, acq, n_code_per);


%% End Matters

% save data with the filter order and the run name in the file
save(['output_run_' run '_filter_' num2str(filter_order) '.mat']);
toc



