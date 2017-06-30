function C=getCsim(cndts,condition_var,A,B,AN, CFs)
C = struct(); % necessary to create structure in parfor
C.cF_i = CFs(cndts(condition_var,1)); %#ok<*PFBNS>
C.sentence_i = A.sentences(cndts(condition_var,2));
C.level= A.level(cndts(condition_var,3));
C.noise_i = B.noiseTypes{cndts(condition_var,4)}; %'SSN';%cs(c_i,4);
C.snr_i = B.SNR(cndts(condition_var,5));
C.ftype_i = AN.fiberType.sr(cndts(condition_var,6));

% resultPostfix = sprintf(resultTxt,       C.cF_i/1e3,  C.sentence_i,  C.noise_i, C.level,    C.snr_i);