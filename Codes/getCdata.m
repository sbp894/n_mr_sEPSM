function C=getCdata(spike_data,condition_var)

paramsIN=spike_data(condition_var).paramsIN;

C = struct(); % necessary to create structure in parfor

C.cF_i=spike_data(condition_var).CF;
C.sentence_i = 1;
C.snr_i = spike_data(condition_var).SNR;
C.noise_i = 'SSN';

C.level=paramsIN.level;