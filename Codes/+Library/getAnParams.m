function paramsIN = getAnParams(duration_A, duration_B, AN,cF_i)

% specify params to be used
paramsIN.durA_msec=duration_A*1E3;
paramsIN.durB_msec=duration_B*1E3;
paramsIN.CF_A_Hz=AN.CF(cF_i);
paramsIN.CF_B_Hz=AN.CF(cF_i);
% Need to include CF_A, CF_B for more generality
paramsIN.MAXspikes=2500;
paramsIN.PSD_LHfreqs_Hz = [1 1; 2 3; 2 6; 4 12; 8 24; 16 48];
paramsIN.UserDelays_usec=[1000 2000];  % Extra delays at which to compute correlations (0 and CD are defaults)
paramsIN.MAXdelay_sec=1;
paramsIN.DELAYbinwidth_sec = 1/50000;
% paramsIN.PSD_LHfreqs_Hz=[0 64; 0 50];  %additional freq ranges to compute CCCenv for
% Can specify which SCpeak and CCCenv to use
% (DEFAULTS: SCpeak='adj: 0-CF';CCCenv='10-300, subBIAS')
% paramsIN.SCpeak_TOUSE='IFFTraw';
% paramsIN.CCCenv_TOUSE='10-300, subBIAS';
end