function [SACSCCs,AvgRate_sps]=SACSCCanal(ST_A_plus,ST_A_minus,ST_B_plus,ST_B_minus,ST_C_plus,ST_C_minus,paramsOUT, winCorr0Add1Mul)

MAXdelay_ind=round(paramsOUT.MAXdelay_sec/paramsOUT.DELAYbinwidth_sec);  % old XLIM=250
% if isfield(ExpControlParams, 'winCorr0Add1Mul')
%     winCorr0Add1Mul=ExpControlParams.winCorr0Add1Mul;
% else
%     winCorr0Add1Mul=1; % 0 for additive, 1 for multiplicative
% end

%% COLUMN 1: CONDITION=A
%% Compute SAC (A+ and A-) %%%%%%%%%%%%%%%%%%%%
% [SAC_A_plus,~,AvgRate_sps(1,1),~] = Library.ShufAutoCorr(ST_A_plus,paramsOUT.DELAYbinwidth_sec,paramsOUT.SCCdur_sec);
% [SAC_A_minus,delays_usec,AvgRate_sps(1,2),~] = Library.ShufAutoCorr(ST_A_minus,paramsOUT.DELAYbinwidth_sec,paramsOUT.SCCdur_sec);

[SAC_A_plus,~,AvgRate_sps(1,1),~] = Library.SAChalf_m(ST_A_plus,paramsOUT.DELAYbinwidth_sec,paramsOUT.SCCdur_sec);
[SAC_A_minus,delays,AvgRate_sps(1,2),~] = Library.SAChalf_m(ST_A_minus,paramsOUT.DELAYbinwidth_sec,paramsOUT.SCCdur_sec);
delays_usec=delays*1e6;

% SAME for all conditions
ZEROind=find(delays_usec==0);
SACinds=(ZEROind-MAXdelay_ind:ZEROind+MAXdelay_ind);
% SACdelays_usec=delays_usec(SACinds);
% SAC limited-data Window CORRECTION
TEMP=linspace(0,1,ZEROind);
if winCorr0Add1Mul
    TEMP1=1./(1-TEMP);
    WindowCORRECTION=[fliplr(TEMP1(2:end)) TEMP1];
    WindowCORRECTION(isinf(WindowCORRECTION))=1/eps;
else
    WindowCORRECTION=[fliplr(TEMP(2:end)) TEMP];
end

% % USED TO TAKE OUT WINDOW CORRECTION to test
% beep
% disp('NO WINDOW CORRECTION APPLIED')
% WindowCORRECTION=zeros(size(WindowCORRECTION));

if winCorr0Add1Mul
    SAC_A_plus=SAC_A_plus.*WindowCORRECTION;
    SAC_A_minus=SAC_A_minus.*WindowCORRECTION;
else
    SAC_A_plus=SAC_A_plus+WindowCORRECTION;
    SAC_A_minus=SAC_A_minus+WindowCORRECTION;
end

% smooth and AVG correlograms
SAC_A_plus=Library.trifilt(SAC_A_plus(SACinds),paramsOUT.TriFiltWidthSAC);
SAC_A_minus=Library.trifilt(SAC_A_minus(SACinds),paramsOUT.TriFiltWidthSAC);
SAC_A_avg=(SAC_A_plus+SAC_A_minus)/2;

%% Compute XpAC (A+/A- and A-/A+) %%%%%%%%%%%%%%%%%%%%
[XpAC_A_plus,~,~,~] = Library.SCCfull_m({ST_A_plus ST_A_minus},paramsOUT.DELAYbinwidth_sec,paramsOUT.SCCdur_sec);
% SAC limited-data Window CORRECTION
if winCorr0Add1Mul
    XpAC_A_plus=XpAC_A_plus.*WindowCORRECTION;
else
    XpAC_A_plus=XpAC_A_plus+WindowCORRECTION;
end
% smooth and AVG correlograms
XpAC_A_plus=Library.trifilt(XpAC_A_plus(SACinds),paramsOUT.TriFiltWidthSAC);
XpAC_A_minus=fliplr(XpAC_A_plus);
XpAC_A_avg=(XpAC_A_plus+XpAC_A_minus)/2;

%% Compute DIFCOR and SUMCOR for Condition A
% % % % DIFCOR_A=Library.trifilt(SAC_A_avg-XpAC_A_avg,paramsOUT.TriFiltWidthDC);  % Diffcorr =  Avg(SAC) - Avg(XpAC)
SUMCOR_A=Library.trifilt((SAC_A_avg+XpAC_A_avg)/2,paramsOUT.TriFiltWidthSC); %Sumcorr = 1/2(Avg(SAC) + Avg(XpAC))

%% Remove TFS artifact (centered at 2*CF) from SUMCOR
% zero out above CF
% FFTtemp=fft((SUMCOR_A-1),paramsOUT.Nfft_psd);
FFTtemp=fft((SUMCOR_A-mean(SUMCOR_A)),paramsOUT.Nfft_psd)/length(SUMCOR_A);
Fs_PSD=1/paramsOUT.DELAYbinwidth_sec;  % Sampling Rate for PSDs
freqVEC=(0:paramsOUT.Nfft_psd-1)/paramsOUT.Nfft_psd*Fs_PSD;
[~,CF_index]=min(abs(freqVEC-paramsOUT.SACSCC_CF_Hz)); % use SACSCC_CF_Hz
FFTadj=zeros(size(FFTtemp));
FFTadj(1:CF_index)=FFTtemp(1:CF_index);
FFTadj((length(FFTtemp)-CF_index+1):end)=FFTtemp((length(FFTtemp)-CF_index+1):end); %keep negative freqs
% % % % adjSC=ifft(FFTadj)+1;
% % % % SUMCORadj_A=real(adjSC(1:length(SUMCOR_A)));
% % % %
% % % % %%%%%SUMCOR IFFT 0-5 Hz%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % % [~,CF_index_5]=min(abs(freqVEC-5)); % use CF_A
% % % % FFTadj_5=zeros(size(FFTtemp));
% % % % FFTadj_5(1:CF_index_5)=FFTtemp(1:CF_index_5);
% % % % FFTadj_5((length(FFTtemp)-CF_index_5+1):end)=FFTtemp((length(FFTtemp)-CF_index_5+1):end); %keep negative freqs
% % % % adjSC_5=ifft(FFTadj_5)+1;
% % % % SUMCORadj_A_5=real(adjSC_5(1:length(SUMCOR_A)));
% % % % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % %
% % % % %%%%%SUMCOR IFFT 0-64 Hz%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % % [~,CF_index_64]=min(abs(freqVEC-64)); % use CF_A
% % % % FFTadj_64=zeros(size(FFTtemp));
% % % % FFTadj_64(1:CF_index_64)=FFTtemp(1:CF_index_64);
% % % % FFTadj_64((length(FFTtemp)-CF_index_64+1):end)=FFTtemp((length(FFTtemp)-CF_index_64+1):end); %keep negative freqs
% % % % adjSC_64=ifft(FFTadj_64)+1;
% % % % SUMCORadj_A_64=real(adjSC_64(1:length(SUMCOR_A)));
% % % % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % %
% % % % %%%%%SUMCOR IFFT 0-300 Hz%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % % [~,CF_index_300]=min(abs(freqVEC-300)); % use CF_A
% % % % FFTadj_300=zeros(size(FFTtemp));
% % % % FFTadj_300(1:CF_index_300)=FFTtemp(1:CF_index_300);
% % % % FFTadj_300((length(FFTtemp)-CF_index_300+1):end)=FFTtemp((length(FFTtemp)-CF_index_300+1):end); %keep negative freqs
% % % % adjSC_300=ifft(FFTadj_300)+1;
% % % % SUMCORadj_A_300=real(adjSC_300(1:length(SUMCOR_A)));
% % % % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % %
% % % % %%%%%SUMCOR IFFT 5-64 Hz%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % % FFTadj_5_64=zeros(size(FFTtemp));
% % % % FFTadj_5_64(CF_index_5:CF_index_64)=FFTtemp(CF_index_5:CF_index_64);
% % % % FFTadj_5_64((length(FFTtemp)-CF_index_64+1):(length(FFTtemp)-CF_index_5+1))=FFTtemp((length(FFTtemp)-CF_index_64+1):(length(FFTtemp)-CF_index_5+1)); %keep negative freqs
% % % % adjSC_5_64=ifft(FFTadj_5_64)+1;
% % % % SUMCORadj_A_5_64=real(adjSC_5_64(1:length(SUMCOR_A)));
% % % % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % %
% % % % %%%%%SUMCOR IFFT 64-300 Hz%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % % FFTadj_64_300=zeros(size(FFTtemp));
% % % % FFTadj_64_300(CF_index_64:CF_index_300)=FFTtemp(CF_index_64:CF_index_300);
% % % % FFTadj_64_300((length(FFTtemp)-CF_index_300+1):(length(FFTtemp)-CF_index_64+1))=FFTtemp((length(FFTtemp)-CF_index_300+1):(length(FFTtemp)-CF_index_64+1)); %keep negative freqs
% % % % adjSC_64_300=ifft(FFTadj_64_300)+1;
% % % % SUMCORadj_A_64_300=real(adjSC_64_300(1:length(SUMCOR_A)));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Store ENVELOPE POWER SPECTRAL DENSITY for Condition A
PSDenv_A=abs(FFTadj);  % use freqVEC to store


%% COLUMN 2: CONDITION=B
% % % % %% Compute SAC (B+ and B-) %%%%%%%%%%%%%%%%%%%%
[SAC_B_plus,~,AvgRate_sps(2,1),~] = Library.SAChalf_m(ST_B_plus,paramsOUT.DELAYbinwidth_sec,paramsOUT.SCCdur_sec);
[SAC_B_minus,~,AvgRate_sps(2,2),~] = Library.SAChalf_m(ST_B_minus,paramsOUT.DELAYbinwidth_sec,paramsOUT.SCCdur_sec);
% SAC limited-data Window CORRECTION
if winCorr0Add1Mul
    SAC_B_plus=SAC_B_plus.*WindowCORRECTION;
    SAC_B_minus=SAC_B_minus.*WindowCORRECTION;
else
    SAC_B_plus=SAC_B_plus+WindowCORRECTION;
    SAC_B_minus=SAC_B_minus+WindowCORRECTION;
end
% smooth and AVG correlograms

SAC_B_plus=Library.trifilt(SAC_B_plus(SACinds),paramsOUT.TriFiltWidthSAC);
SAC_B_minus=Library.trifilt(SAC_B_minus(SACinds),paramsOUT.TriFiltWidthSAC);
SAC_B_avg=(SAC_B_plus+SAC_B_minus)/2;

%% Compute XpAC (B+/B- and B-/B+) %%%%%%%%%%%%%%%%%%%%
[XpAC_B_plus,~,~,~] = Library.SCCfull_m({ST_B_plus ST_B_minus},paramsOUT.DELAYbinwidth_sec,paramsOUT.SCCdur_sec);
% SAC limited-data Window CORRECTION
if winCorr0Add1Mul
    XpAC_B_plus=XpAC_B_plus.*WindowCORRECTION;
else
    XpAC_B_plus=XpAC_B_plus+WindowCORRECTION;
end
% smooth and AVG correlograms
XpAC_B_plus=Library.trifilt(XpAC_B_plus(SACinds),paramsOUT.TriFiltWidthSAC);
XpAC_B_minus=fliplr(XpAC_B_plus);
XpAC_B_avg=(XpAC_B_plus+XpAC_B_minus)/2;

%% Compute DIFCOR and SUMCOR for Condition B
% % % % DIFCOR_B=Library.trifilt(SAC_B_avg-XpAC_B_avg,paramsOUT.TriFiltWidthDC);  % Diffcorr =  Avg(SAC) - Avg(XpAC)
SUMCOR_B=Library.trifilt((SAC_B_avg+XpAC_B_avg)/2,paramsOUT.TriFiltWidthSC); %Sumcorr = 1/2(Avg(SAC) + Avg(XpAC))

%% Remove TFS artifact (centered at 2*CF) from SUMCOR
% zero out above CF
% FFTtemp=fft((SUMCOR_B-1),paramsOUT.Nfft_psd);
FFTtemp=fft((SUMCOR_B-mean(SUMCOR_B)),paramsOUT.Nfft_psd)/length(SUMCOR_B);
[~,CF_index]=min(abs(freqVEC-paramsOUT.SACSCC_CF_Hz)); % use SACSCC_CF_Hz
FFTadj=zeros(size(FFTtemp));
FFTadj(1:CF_index)=FFTtemp(1:CF_index);
FFTadj((length(FFTtemp)-CF_index+1):end)=FFTtemp((length(FFTtemp)-CF_index+1):end); %keep negative freqs
% % % % adjSC=ifft(FFTadj)+1;
% % % % SUMCORadj_B=real(adjSC(1:length(SUMCOR_B)));
% % % %
% % % % %%%%%SUMCOR IFFT 0-5 Hz%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % % [~,CF_index_5]=min(abs(freqVEC-5)); % use CF_B
% % % % FFTadj_5=zeros(size(FFTtemp));
% % % % FFTadj_5(1:CF_index_5)=FFTtemp(1:CF_index_5);
% % % % FFTadj_5((length(FFTtemp)-CF_index_5+1):end)=FFTtemp((length(FFTtemp)-CF_index_5+1):end); %keep negative freqs
% % % % adjSC_5=ifft(FFTadj_5)+1;
% % % % SUMCORadj_B_5=real(adjSC_5(1:length(SUMCOR_B)));
% % % % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % %
% % % % [~,CF_index_64]=min(abs(freqVEC-64)); % use CF_A
% % % % FFTadj_64=zeros(size(FFTtemp));
% % % % FFTadj_64(1:CF_index_64)=FFTtemp(1:CF_index_64);
% % % % FFTadj_64((length(FFTtemp)-CF_index_64+1):end)=FFTtemp((length(FFTtemp)-CF_index_64+1):end); %keep negative freqs
% % % % adjSC_64=ifft(FFTadj_64)+1;
% % % % SUMCORadj_B_64=real(adjSC_64(1:length(SUMCOR_B)));
% % % %
% % % % [~,CF_index_300]=min(abs(freqVEC-300)); % use CF_A
% % % % FFTadj_300=zeros(size(FFTtemp));
% % % % FFTadj_300(1:CF_index_300)=FFTtemp(1:CF_index_300);
% % % % FFTadj_300((length(FFTtemp)-CF_index_300+1):end)=FFTtemp((length(FFTtemp)-CF_index_300+1):end); %keep negative freqs
% % % % adjSC_300=ifft(FFTadj_300)+1;
% % % % SUMCORadj_B_300=real(adjSC_300(1:length(SUMCOR_B)));
% % % %
% % % % %%%%%SUMCOR IFFT 5-64 Hz%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % % FFTadj_5_64=zeros(size(FFTtemp));
% % % % FFTadj_5_64(CF_index_5:CF_index_64)=FFTtemp(CF_index_5:CF_index_64);
% % % % FFTadj_5_64((length(FFTtemp)-CF_index_64+1):(length(FFTtemp)-CF_index_5+1))=FFTtemp((length(FFTtemp)-CF_index_64+1):(length(FFTtemp)-CF_index_5+1)); %keep negative freqs
% % % % adjSC_5_64=ifft(FFTadj_5_64)+1;
% % % % SUMCORadj_B_5_64=real(adjSC_5_64(1:length(SUMCOR_B)));
% % % % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % %
% % % % FFTadj_64_300=zeros(size(FFTtemp));
% % % % FFTadj_64_300(CF_index_64:CF_index_300)=FFTtemp(CF_index_64:CF_index_300);
% % % % FFTadj_64_300((length(FFTtemp)-CF_index_300+1):(length(FFTtemp)-CF_index_64+1))=FFTtemp((length(FFTtemp)-CF_index_300+1):(length(FFTtemp)-CF_index_64+1)); %keep negative freqs
% % % % adjSC_64_300=ifft(FFTadj_64_300)+1;
% % % % SUMCORadj_B_64_300=real(adjSC_64_300(1:length(SUMCOR_B)));

%% Compute ENVELOPE POWER SPECTRAL DENSITY for Condition B
PSDenv_B=abs(FFTadj);  % use freqVEC to store
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% COLUMN 3: CONDITION=C
%% Compute SAC (C+ and C-) %%%%%%%%%%%%%%%%%%%%
[SAC_C_plus,~,AvgRate_sps(3,1),~] = Library.SAChalf_m(ST_C_plus,paramsOUT.DELAYbinwidth_sec,paramsOUT.SCCdur_sec);
[SAC_C_minus,~,AvgRate_sps(3,2),~] = Library.SAChalf_m(ST_C_minus,paramsOUT.DELAYbinwidth_sec,paramsOUT.SCCdur_sec);
% SAC limited-data Window CORRECTION
if winCorr0Add1Mul
    SAC_C_plus=SAC_C_plus.*WindowCORRECTION;
    SAC_C_minus=SAC_C_minus.*WindowCORRECTION;
else
    SAC_C_plus=SAC_C_plus+WindowCORRECTION;
    SAC_C_minus=SAC_C_minus+WindowCORRECTION;
end
% smooth and AVG correlograms
SAC_C_plus=Library.trifilt(SAC_C_plus(SACinds),paramsOUT.TriFiltWidthSAC);
SAC_C_minus=Library.trifilt(SAC_C_minus(SACinds),paramsOUT.TriFiltWidthSAC);
SAC_C_avg=(SAC_C_plus+SAC_C_minus)/2;

%% Compute XpAC (C+/C- and C-/C+) %%%%%%%%%%%%%%%%%%%%
[XpAC_C_plus,~,~,~] = Library.SCCfull_m({ST_C_plus ST_C_minus},paramsOUT.DELAYbinwidth_sec,paramsOUT.SCCdur_sec);
% SAC limited-data Window CORRECTION
if winCorr0Add1Mul
    XpAC_C_plus=XpAC_C_plus.*WindowCORRECTION;
else
    XpAC_C_plus=XpAC_C_plus+WindowCORRECTION;
end
% smooth and AVG correlograms
XpAC_C_plus=Library.trifilt(XpAC_C_plus(SACinds),paramsOUT.TriFiltWidthSAC);
XpAC_C_minus=fliplr(XpAC_C_plus);
XpAC_C_avg=(XpAC_C_plus+XpAC_C_minus)/2;

%% Compute DIFCOR and SUMCOR for Condition C
% % % % DIFCOR_C=Library.trifilt(SAC_C_avg-XpAC_C_avg,paramsOUT.TriFiltWidthDC);  % Diffcorr =  Avg(SAC) - Avg(XpAC)
SUMCOR_C=Library.trifilt((SAC_C_avg+XpAC_C_avg)/2,paramsOUT.TriFiltWidthSC); %Sumcorr = 1/2(Avg(SAC) + Avg(XpAC))

%% Remove TFS artifact (centered at 2*CF) from SUMCOR
% zero out above CF
% FFTtemp=fft((SUMCOR_C-1),paramsOUT.Nfft_psd);
FFTtemp=fft((SUMCOR_C-mean(SUMCOR_C)),paramsOUT.Nfft_psd)/length(SUMCOR_C);
[~,CF_index]=min(abs(freqVEC-paramsOUT.SACSCC_CF_Hz)); % use SACSCC_CF_Hz
FFTadj=zeros(size(FFTtemp));
FFTadj(1:CF_index)=FFTtemp(1:CF_index);
FFTadj((length(FFTtemp)-CF_index+1):end)=FFTtemp((length(FFTtemp)-CF_index+1):end); %keep negative freqs
% % % % adjSC=ifft(FFTadj)+1;
% % % % SUMCORadj_C=real(adjSC(1:length(SUMCOR_C)));
% % % %
% % % % %%%%%SUMCOR IFFT 0-5 Hz%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % % [~,CF_index_5]=min(abs(freqVEC-5)); % use CF_C
% % % % FFTadj_5=zeros(size(FFTtemp));
% % % % FFTadj_5(1:CF_index_5)=FFTtemp(1:CF_index_5);
% % % % FFTadj_5((length(FFTtemp)-CF_index_5+1):end)=FFTtemp((length(FFTtemp)-CF_index_5+1):end); %keep negative freqs
% % % % adjSC_5=ifft(FFTadj_5)+1;
% % % % SUMCORadj_C_5=real(adjSC_5(1:length(SUMCOR_C)));
% % % % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % %
% % % % [~,CF_index_64]=min(abs(freqVEC-64)); % use CF_A
% % % % FFTadj_64=zeros(size(FFTtemp));
% % % % FFTadj_64(1:CF_index_64)=FFTtemp(1:CF_index_64);
% % % % FFTadj_64((length(FFTtemp)-CF_index_64+1):end)=FFTtemp((length(FFTtemp)-CF_index_64+1):end); %keep negative freqs
% % % % adjSC_64=ifft(FFTadj_64)+1;
% % % % SUMCORadj_C_64=real(adjSC_64(1:length(SUMCOR_C)));
% % % %
% % % % [~,CF_index_300]=min(abs(freqVEC-300)); % use CF_A
% % % % FFTadj_300=zeros(size(FFTtemp));
% % % % FFTadj_300(1:CF_index_300)=FFTtemp(1:CF_index_300);
% % % % FFTadj_300((length(FFTtemp)-CF_index_300+1):end)=FFTtemp((length(FFTtemp)-CF_index_300+1):end); %keep negative freqs
% % % % adjSC_300=ifft(FFTadj_300)+1;
% % % % SUMCORadj_C_300=real(adjSC_300(1:length(SUMCOR_C)));
% % % %
% % % % %%%%%SUMCOR IFFT 5-64 Hz%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % % FFTadj_5_64=zeros(size(FFTtemp));
% % % % FFTadj_5_64(CF_index_5:CF_index_64)=FFTtemp(CF_index_5:CF_index_64);
% % % % FFTadj_5_64((length(FFTtemp)-CF_index_64+1):(length(FFTtemp)-CF_index_5+1))=FFTtemp((length(FFTtemp)-CF_index_64+1):(length(FFTtemp)-CF_index_5+1)); %keep negative freqs
% % % % adjSC_5_64=ifft(FFTadj_5_64)+1;
% % % % SUMCORadj_C_5_64=real(adjSC_5_64(1:length(SUMCOR_C)));
% % % % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % %
% % % % FFTadj_64_300=zeros(size(FFTtemp));
% % % % FFTadj_64_300(CF_index_64:CF_index_300)=FFTtemp(CF_index_64:CF_index_300);
% % % % FFTadj_64_300((length(FFTtemp)-CF_index_300+1):(length(FFTtemp)-CF_index_64+1))=FFTtemp((length(FFTtemp)-CF_index_300+1):(length(FFTtemp)-CF_index_64+1)); %keep negative freqs
% % % % adjSC_64_300=ifft(FFTadj_64_300)+1;
% % % % SUMCORadj_C_64_300=real(adjSC_64_300(1:length(SUMCOR_C)));

%% Compute ENVELOPE POWER SPECTRAL DENSITY for Condition C
PSDenv_C=abs(FFTadj);  % use freqVEC to store
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% COLUMN 4: CONDITION=AC (clean speech, noisy speech)
%% Compute SCC (A+/C+ and A-/C-) %%%%%%%%%%%%%%%%%%%%
% % % % [SCC_AC_plus,~,~,~] = Library.ShufCrossCorr({ST_A_plus ST_C_plus},paramsOUT.DELAYbinwidth_sec,paramsOUT.SCCdur_sec);
% % % % [SCC_AC_minus,~,~,~] = Library.ShufCrossCorr({ST_A_minus ST_C_minus},paramsOUT.DELAYbinwidth_sec,paramsOUT.SCCdur_sec);
% % % % % SAC limited-data Window CORRECTION
% % % % SCC_AC_plus=SCC_AC_plus.*WindowCORRECTION;
% % % % SCC_AC_minus=SCC_AC_minus.*WindowCORRECTION;
% % % % % SCC_AC_plus=SCC_AC_plus+WindowCORRECTION;
% % % % % SCC_AC_minus=SCC_AC_minus+WindowCORRECTION;
% % % % % smooth and AVG correlograms
% % % % SCC_AC_plus=Library.trifilt(SCC_AC_plus(SACinds),paramsOUT.TriFiltWidthSAC);
% % % % SCC_AC_minus=Library.trifilt(SCC_AC_minus(SACinds),paramsOUT.TriFiltWidthSAC);
% % % % SCC_AC_avg=(SCC_AC_plus+SCC_AC_minus)/2; % That will be plotted in first row of plot
% % % %
% % % % %% Compute XpCC (A+/B- and A-/B+) %%%%%%%%%%%%%%%%%%%%
% % % % [XpCC_AC_plus,~,~,~] = Library.ShufCrossCorr({ST_A_plus ST_B_minus},paramsOUT.DELAYbinwidth_sec,paramsOUT.SCCdur_sec);
% % % % [XpCC_AC_minus,~,~,~] = Library.ShufCrossCorr({ST_A_minus ST_B_plus},paramsOUT.DELAYbinwidth_sec,paramsOUT.SCCdur_sec);
% % % % % SAC limited-data Window CORRECTION
% % % % XpCC_AC_plus=XpCC_AC_plus.*WindowCORRECTION;
% % % % XpCC_AC_minus=XpCC_AC_minus.*WindowCORRECTION;
% % % % % XpCC_AC_plus=XpCC_AC_plus+WindowCORRECTION;
% % % % % XpCC_AC_minus=XpCC_AC_minus+WindowCORRECTION;
% % % % % smooth and AVG correlograms
% % % % XpCC_AC_plus=Library.trifilt(XpCC_AC_plus(SACinds),paramsOUT.TriFiltWidthSAC);
% % % % XpCC_AC_minus=Library.trifilt(XpCC_AC_minus(SACinds),paramsOUT.TriFiltWidthSAC);
% % % % XpCC_AC_avg=(XpCC_AC_plus+XpCC_AC_minus)/2; % that will be plotted in second row of plot
% % % %
% % % % %% Compute DIFCOR and SUMCOR for Condition AC
% % % % DIFCOR_AC=Library.trifilt(SCC_AC_avg-XpCC_AC_avg,paramsOUT.TriFiltWidthDC);  % Diffcorr =  Avg(SCC) - Avg(XpCC) that will be plotted in third row of plot
% % % % SUMCOR_AC=Library.trifilt((SCC_AC_avg+XpCC_AC_avg)/2,paramsOUT.TriFiltWidthSC); %Sumcorr = 1/2(Avg(SCC) + Avg(XpCC))
% % % %
% % % % %% Remove TFS artifact (centered at 2*CF) from SUMCOR
% % % % % zero out above CF
% % % % % FFTtemp=fft((SUMCOR_AC-1),paramsOUT.Nfft_psd);
% % % % FFTtemp=fft((SUMCOR_AC-mean(SUMCOR_AC)),paramsOUT.Nfft_psd)/length(SUMCOR_AC);
% % % % [~,CF_index]=min(abs(freqVEC-paramsOUT.SACSCC_CF_Hz)); % use min([CF_A CF_B])
% % % % FFTadj=zeros(size(FFTtemp));
% % % % FFTadj(1:CF_index)=FFTtemp(1:CF_index);
% % % % FFTadj((length(FFTtemp)-CF_index+1):end)=FFTtemp((length(FFTtemp)-CF_index+1):end); %keep negative freqs
% % % % adjSC=ifft(FFTadj)+1;
% % % % SUMCORadj_AC=real(adjSC(1:length(SUMCOR_AC))); % will be plotted in
% % % %
% % % % %%%%%SUMCOR IFFT 0-5 Hz%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % % [~,CF_index_5]=min(abs(freqVEC-5)); % use CF_B
% % % % FFTadj_5=zeros(size(FFTtemp));
% % % % FFTadj_5(1:CF_index_5)=FFTtemp(1:CF_index_5);
% % % % FFTadj_5((length(FFTtemp)-CF_index_5+1):end)=FFTtemp((length(FFTtemp)-CF_index_5+1):end); %keep negative freqs
% % % % adjSC_5=ifft(FFTadj_5)+1;
% % % % SUMCORadj_AC_5=real(adjSC_5(1:length(SUMCOR_AC)));
% % % % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % %
% % % % [~,CF_index_64]=min(abs(freqVEC-64)); % use CF_A
% % % % FFTadj_64=zeros(size(FFTtemp));
% % % % FFTadj_64(1:CF_index_64)=FFTtemp(1:CF_index_64);
% % % % FFTadj_64((length(FFTtemp)-CF_index_64+1):end)=FFTtemp((length(FFTtemp)-CF_index_64+1):end); %keep negative freqs
% % % % adjSC_64=ifft(FFTadj_64)+1;
% % % % SUMCORadj_AC_64=real(adjSC_64(1:length(SUMCOR_AC)));
% % % %
% % % % [~,CF_index_300]=min(abs(freqVEC-300)); % use CF_A
% % % % FFTadj_300=zeros(size(FFTtemp));
% % % % FFTadj_300(1:CF_index_300)=FFTtemp(1:CF_index_300);
% % % % FFTadj_300((length(FFTtemp)-CF_index_300+1):end)=FFTtemp((length(FFTtemp)-CF_index_300+1):end); %keep negative freqs
% % % % adjSC_300=ifft(FFTadj_300)+1;
% % % % SUMCORadj_AC_300=real(adjSC_300(1:length(SUMCOR_AC)));
% % % %
% % % % %%%%%SUMCOR IFFT 5-64 Hz%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % % FFTadj_5_64=zeros(size(FFTtemp));
% % % % FFTadj_5_64(CF_index_5:CF_index_64)=FFTtemp(CF_index_5:CF_index_64);
% % % % FFTadj_5_64((length(FFTtemp)-CF_index_64+1):(length(FFTtemp)-CF_index_5+1))=FFTtemp((length(FFTtemp)-CF_index_64+1):(length(FFTtemp)-CF_index_5+1)); %keep negative freqs
% % % % adjSC_5_64=ifft(FFTadj_5_64)+1;
% % % % SUMCORadj_AC_5_64=real(adjSC_5_64(1:length(SUMCOR_AC)));
% % % % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % %
% % % % FFTadj_64_300=zeros(size(FFTtemp));
% % % % FFTadj_64_300(CF_index_64:CF_index_300)=FFTtemp(CF_index_64:CF_index_300);
% % % % FFTadj_64_300((length(FFTtemp)-CF_index_300+1):(length(FFTtemp)-CF_index_64+1))=FFTtemp((length(FFTtemp)-CF_index_300+1):(length(FFTtemp)-CF_index_64+1)); %keep negative freqs
% % % % adjSC_64_300=ifft(FFTadj_64_300)+1;
% % % % SUMCORadj_AC_64_300=real(adjSC_64_300(1:length(SUMCOR_AC)));
% % % %
% % % %
% % % % %% Compute ENVELOPE CROSS SPECTRAL DENSITY for Condition AC
% % % % PSDenv_AC=abs(FFTtemp);
% % % % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % %
% % % % %% Compute Spectral Density of DIFCOR_A
% % % % FFTtemp=fft((DIFCOR_A),paramsOUT.Nfft_psd);
% % % % PSDtfs_A=abs(FFTtemp);
% % % %
% % % % %% Compute Spectral Density of DIFCOR_B
% % % % FFTtemp=fft((DIFCOR_B),paramsOUT.Nfft_psd);
% % % % PSDtfs_B=abs(FFTtemp);
% % % %
% % % % %% Compute Spectral Density of DIFCOR_C
% % % % FFTtemp=fft((DIFCOR_C),paramsOUT.Nfft_psd);
% % % % PSDtfs_C=abs(FFTtemp);
% % % %
% % % % %% Compute Spectral Density of DIFCOR_AC
% % % % FFTtemp=fft((DIFCOR_AC),paramsOUT.Nfft_psd);
% % % % PSDtfs_AC=abs(FFTtemp);

%% Store all params, output for return;
% SACSCCs=struct('delays_usec',SACdelays_usec,'SAC_A_avg',SAC_A_avg,'SAC_B_avg',SAC_B_avg,'SAC_C_avg',SAC_C_avg,'SCC_AC_avg',SCC_AC_avg,'XpAC_A_avg',XpAC_A_avg, ...
%     'XpAC_B_avg',XpAC_B_avg,'XpAC_C_avg',XpAC_C_avg,'XpCC_AC_avg',XpCC_AC_avg,'SUMCOR_A',SUMCOR_A, ...
%     'SUMCOR_B',SUMCOR_B,'SUMCOR_C',SUMCOR_C,'SUMCOR_AC',SUMCOR_AC,'SUMCORadj_A',SUMCORadj_A,'SUMCORadj_B',SUMCORadj_B,'SUMCORadj_C',SUMCORadj_C,'SUMCORadj_AC',SUMCORadj_AC, ...
%     'DIFCOR_A',DIFCOR_A,'DIFCOR_B',DIFCOR_B,'DIFCOR_C',DIFCOR_C,'DIFCOR_AC',DIFCOR_AC, ...
%     'SUMCORadj_A_64',SUMCORadj_A_64,'SUMCORadj_B_64',SUMCORadj_B_64,'SUMCORadj_C_64',SUMCORadj_C_64,'SUMCORadj_AC_64',SUMCORadj_AC_64, ...
%     'SUMCORadj_A_5',SUMCORadj_A_5,'SUMCORadj_B_5',SUMCORadj_B_5,'SUMCORadj_C_5',SUMCORadj_C_5,'SUMCORadj_AC_5',SUMCORadj_AC_5, ...
%     'SUMCORadj_A_300',SUMCORadj_A_300,'SUMCORadj_B_300',SUMCORadj_B_300,'SUMCORadj_C_300',SUMCORadj_C_300,'SUMCORadj_AC_300',SUMCORadj_AC_300, ...
%     'SUMCORadj_A_64_300',SUMCORadj_A_64_300,'SUMCORadj_B_64_300',SUMCORadj_B_64_300,'SUMCORadj_C_64_300',SUMCORadj_C_64_300,'SUMCORadj_AC_64_300',SUMCORadj_AC_64_300, ...
%     'SUMCORadj_A_5_64',SUMCORadj_A_5_64,'SUMCORadj_B_5_64',SUMCORadj_B_5_64,'SUMCORadj_C_5_64',SUMCORadj_C_5_64,'SUMCORadj_AC_5_64',SUMCORadj_AC_5_64,'PSDenv_A',PSDenv_A,...
%     'PSDenv_B',PSDenv_B,'PSDenv_C',PSDenv_C,'PSDenv_AC',PSDenv_AC,'PSDtfs_A',PSDtfs_A,'PSDtfs_B',PSDtfs_B,'PSDtfs_C',PSDtfs_C,'PSDtfs_AC',PSDtfs_AC,'PSD_freqVEC',freqVEC,...
%     'freqVEC',freqVEC);

% SACSCCs=struct('PSDenv_A',PSDenv_A,'PSDenv_B',PSDenv_B,'PSDenv_C',PSDenv_C,'PSD_freqVEC',freqVEC);

SACSCCs=struct('delays_usec',delays_usec(SACinds), 'SUMCOR_A',SUMCOR_A,'SUMCOR_B',SUMCOR_B,'SUMCOR_C',SUMCOR_C,'PSDenv_A',PSDenv_A,'PSDenv_B',PSDenv_B,'PSDenv_C',PSDenv_C,'PSD_freqVEC',freqVEC);

return;

