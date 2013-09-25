%% MECH 6970 Lab 2, Part 1
% Robert Cofield, et al.
% 
% due 2013-09-30
% 
% *Run this file in its directory*
%    - that's unrobust, I know, but it's late...
%
% Prerequisites
%   - mydate package (fileexchange)
%   - goGPS (gogps-project.org) - I did add with subfolders
%   - rgc_matlab (github.com/ozymandium/rgc_matlab)
%         - just need to add top level to path
% 
genutil.ccc

%% Configurating
% UTC Time datevec from first GPZDA message of second data file (*.ubx)
dtvec = [2013,09,10,17,37,34]; % here's the hard-coding. oh well.

%% Time Figuring Outing
yrstr = num2str(dtvec(1));
doy = timeutil.datevec2doy(dtvec); % day of year, from rgc_matlab
% with this we can get CDDIS ephems, but want precise ephem
% get gps week and seconds into week
[wk,sec] = mydategps(mydatenum(dtvec)); % you'll need the mydate package installed
% day of week (0=Sun, 6=Sat)
dow = weekday(datestr(dtvec))-1;

%% CDDIS Ephemeris Getting & Parsing
% find the internet address for CDDIS ephem
% cddis = ftp('cddis.gsfc.nasa.gov'); % ftp object for the server
% filename
cddis_addr = ['pub/gps/data/daily/', yrstr, '/brdc/'];
cddis_fname = ['brdc' sprintf('%3.3d',doy) '0.' yrstr(3:4) 'n.Z'];
mkdir('../data/tmp');
% the paths here probably won't work on Windows .. Meh.
system(['python ../data/download_cddis_data.py ' cddis_addr cddis_fname ' ../data/tmp/']);
cddis_fname =  cddis_fname(1:end-2);
[ephem,~] = RINEX_get_nav(['../data/tmp/' cddis_fname]);

%% JPL Precise Ephemeris -- Unfinished.
% % find the internet address for JPL precise ephem
% jpl_addr = ['http://igscb.jpl.nasa.gov/igscb/product/', num2str(wk), '/']; % folder
% jpl_addr = [jpl_addr, 'igr', num2str(wk*10+dow), '.sp3.Z']; % filename
% % urlwrite(jpl_addr, 'precise_ephem_jpl.rinex');

%% SV Position Calculation
