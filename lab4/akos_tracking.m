function [trackRes] = akos_tracking(fileid, acq)
% implement akos' tracking algorithm for preliminary data bit decoding

svidx_ = 1;

fL1 = 154*10.23e6; % L1 frequency,  1.5754e+09 Hz
fs = 16.3676e6; % sampling frequency
fIF = 4.1304e6; % intermediate frequency
fD = 50; % Data message chipping frequency, Hz
fCA = 1.023e6; % C/A code frequency, Hz

TL1 = 1/fL1; % L1 carrier period
Tca_chip = 1.023e-6; % L1 C/A chip period
Tca = 1.0e-3; % L1 C/A sequence period
Td_chip = 1/fD;
Ts = 1/fs; % sampling period
Tid = Tca; % integrate & dump period (sec) - use the C/A code period (1ms)
% Tnav_chip = 20e-3; % period of a single nav msg chip
% Tnav_frame = Tnav_chip*1500; % period of the entire nav msg frame
% Tnav_subframe = Tnav_frame/5; % period of a nav msg subframe

stop_time = 60; % how long to go in sec

% preamble = bin2dec('10001011');  

%% Filter calculations
% use Akos' filter
% !!! WE NEED TO DO THIS OURSELVES !!!

% % % DLL 
dll_noise_bw = 2; % (Hz)
dll_damping = 0.7; % ratio
[t1,t2] = calcLoopCoef(dll_noise_bw, dll_damping, 1.0);
dll_filt = struct(...
  'tau1', t1, ...
  'tau2', t2 ...
  );

% % % PLL 
pll_noise_bw = 25; % Hz
pll_damping = 0.7; % ratio
[t1,t2] = calcLoopCoef(pll_noise_bw, pll_damping, 0.25);
pll_filt = struct(...
  'tau1', t1, ...
  'tau2', t2 ...
  );


%% Tracking Loop

n_code_per = stop_time/Tid; % how many data points we'll have in the end
if rem(n_code_per,1), error('Stop Time Needs to be multiple of 0.001'); end

% needed data
facq = acq.fdopp+fIF;
codePhase0 = round( 1023*fs/fCA - acq.tau_samples );

% settings
dll_corr_spc = 0.5; % correlator spacing on DLL is half C/A code chip

% % Initialize outputs of tracking
trackRes.IE = zeros(1,n_code_per);
trackRes.IP = zeros(1,n_code_per);
trackRes.IL = zeros(1,n_code_per);
trackRes.QE = zeros(1,n_code_per);
trackRes.QP = zeros(1,n_code_per);
trackRes.QL = zeros(1,n_code_per);
trackRes.fcode = inf(1,n_code_per);
trackRes.codeErr = inf(1,n_code_per);
trackRes.codeNCO = inf(1,n_code_per);
trackRes.fcarr = inf(1,n_code_per);
trackRes.carrErr = inf(1,n_code_per);
trackRes.carrNCO = inf(1,n_code_per);
repmat(trackRes,1,acq.nsv);

hwb = waitbar(0,'Tracking Initialized'); % waitbar

% Loop over each satellite
for ch = 1:acq.nsv
  
  fseek(fileid, codePhase0(ch)-1 , 'bof'); % return data file to code phase offset
  
  prn = genprn(acq.svs(ch),1023,[-1 1]);
  prn = [ prn(1023) prn prn(1) ]; % make it possible to do early and late
  
  % % initial states
  fcode = fCA; % NCO C/A frequency starts off
  res_codephase = 0; % chips
  fcarr = facq(ch);
  res_carrphase = 0;
  % Don't know what these do yet
  oldCodeNCO = 0;
  oldCodeErr = 0;
  oldCarrNCO = 0;
  oldCarrErr = 0;
  
  % loop over each I&D period
  for k = 1:n_code_per
    
    % update progress bar
    if ~rem(k,50)
      waitbar((k/n_code_per), hwb, ...
        ['Tracking Channel ' num2str(ch) '/' num2str(acq.nsv) ' - PRN' num2str(acq.svs(ch)) ' : ' num2str(k) '/' num2str(n_code_per)]);
    end
    
    % % get new data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    d_codePhase = fcode/fs; % update the phase step using the current code frequency
    % determine how many bytes to read from file
    blksz = ceil((1023-res_codephase)/d_codePhase);
    [signal,nsamples] = fread(fileid, blksz, 'int8');
    if nsamples~=blksz, error('Cannot read correct size'); end
    signal = signal';
    
    % % % Get E, P, L info %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % % get index within early code vector
    tcode = (res_codephase-dll_corr_spc) ...
            : d_codePhase : ...
            ((blksz-1)*d_codePhase + res_codephase - dll_corr_spc);
    tcode = ceil(tcode)+1;
    Ecode = prn(tcode);
    % % get index within Late code vector
    tcode = (res_codephase+dll_corr_spc) ...
            : d_codePhase : ...
            ((blksz-1)*d_codePhase + res_codephase + dll_corr_spc);
    tcode = ceil(tcode)+1;
    Lcode = prn(tcode);
    % % get index within Prompt code vector
    tcode_ = res_codephase ...
            : d_codePhase : ...
            ((blksz-1)*d_codePhase + res_codephase);
    tcode = ceil(tcode_)+1;
    Pcode = prn(tcode);
    % % recalculate code phase residual
    res_codephase = (tcode_(blksz) + d_codePhase) - 1023;
    clear tcode tcode_
    
    % % generate carrier frequency %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    time = (0:blksz)./fs;
    % generate arguments to the sin/cos functions
    arg = ((fcarr*2*pi)*time+res_carrphase);
    % recalculate the carrier phase residual
    res_carrphase = rem(arg(blksz+1),(2*pi));
    % calculate the sin/cos
    sin_carr = sin(arg(1:blksz));
    cos_carr = cos(arg(1:blksz));
    clear arg
    
    % % % Accumulators %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % % In-phase
    I = sin_carr .* signal;
    IE = sum( Ecode .* I );
    IP = sum( Pcode .* I );
    IL = sum( Lcode .* I );
    % store data
    trackRes(ch).IE(k) = IE;
    trackRes(ch).IP(k) = IP;
    trackRes(ch).IL(k) = IL;
    % % Quatrature-phase
    Q = cos_carr .* signal;
    QE = sum( Ecode .* Q );
    QP = sum( Pcode .* Q );
    QL = sum( Lcode .* Q );
    % store data
    trackRes(ch).QE(k) = QE;
    trackRes(ch).QP(k) = QP;
    trackRes(ch).QL(k) = QL;
    
    % % Find PLL Error & update carrier NCO %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    carrErr = atan(QP/IP) / (2*pi);
    % !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    % update carrier NCO with filter
    % !!! This is where we implement filters of other orders !!!!!!!!!!!!!!!!!!!
    carrNCO = oldCarrNCO + ...
              ( pll_filt.tau2 / pll_filt.tau1 ) * (carrErr-oldCarrErr) + ...
              carrErr*(Tid/pll_filt.tau1);
    % !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    % save data for next timestep
    oldCarrNCO = carrNCO;
    oldCarrErr = carrErr;
    % modify carrier frequency based on the carrier NCO
    fcarr = facq(ch) + carrNCO;
    % store to results
    trackRes(ch).fcarr(k) = fcarr;
    trackRes(ch).carrErr(k) = carrErr;
    trackRes(ch).carrNCO(k) = carrNCO;
    
    % % Find DLL Error & update code NCO %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    early = sqrt( IE^2 + QE^2 );
    late  = sqrt( IL^2 + QL^2 );
    codeErr = (early-late) / (early+late);
    % !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    % update code NCO with filter
    % !!! This is where we implement filters of other orders !!!!!!!!!!!!!!!!!!!
    codeNCO = oldCodeNCO + ...
              ( dll_filt.tau2 / dll_filt.tau1 ) * (codeErr-oldCodeErr) + ...
              codeErr*(Tid/dll_filt.tau1);
    % !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    % save data for next timestep
    oldCodeNCO = codeNCO;
    oldCodeErr = codeErr;
    % modify code frequency based on the the code NCO
    fcode = fCA - codeNCO;
    % store to results
    trackRes(ch).fcode(k) = fcode;   
    trackRes(ch).codeErr(k) = codeErr;
    trackRes(ch).codeNCO(k) = codeNCO;
    
  end
  
end

close(hwb);


end