%% IMPORTANT
% Important to note that order of SummaryVar and SpikeStimulusData are not the same
% NEL data files have attn value. I have assumed 120 dB is the level
% corresponding to 0 dB atten. Hence the line
% "cur_attn=120-SummaryVar(ind2updateInSummaryVar(attn_var)).SPL;".

%%
function SummaryVar=get_revcor(SpikeStimulusData, SummaryVar)

% SummaryVar(file_var).AcousticSNR SummaryVar(file_var).CF_revcor
%%

% Use Stim0dB_N_P/ Stim0dB_N_N to get REVCOR for +6 and 0 dB SN
window_dur=.05; % 50 ms
cfsSUMMary=extractfield(SummaryVar,'BF_TC')';
snrSUMMary=extractfield(SummaryVar,'SNR')';
splSUMMARY=extractfield(SummaryVar,'SPL')';
SNR2use=[-6, 0, 6];

spike_data=SpikeStimulusData.spike_data;
StimsFNames=SpikeStimulusData.StimsFNames;
SummaryVar(1).h_revcor=[];

ind2updateInSummaryVar=zeros(length(spike_data),1);
temp=SummaryVar;

for i=1:length(spike_data)
    cur_snr=spike_data(i).SNR;
    cur_cf=spike_data(i).CF;
    cur_spl=spike_data(i).SPL;
    xx=find((cfsSUMMary==cur_cf & snrSUMMary==cur_snr & splSUMMARY==cur_spl)==1);
    if ~isempty(xx)
        ind2updateInSummaryVar(i)=xx;
        temp(i)=SummaryVar(ind2updateInSummaryVar(i));
%         temp(i).pass1fail0=1;
    else
        ind2updateInSummaryVar(i)=nan;
%         temp(i).pass1fail0=0;
    end
end


parfor i=1:length(spike_data)
    if ~isnan(ind2updateInSummaryVar(i))
        cur_snr=spike_data(i).SNR;
        
        %     cur_cf=spike_data(i).CF;
        cur_spl=spike_data(i).SPL;
        %     ind2updateInSummaryVar=find((cfsSUMMary==cur_cf & snrSUMMary==cur_snr & splSUMMARY==cur_spl)==1);
        if ismember(cur_snr,SNR2use)
            
            cur_attn=120-cur_spl;
            
            %         if cur_snr==6
            %             check0dBInd=find((cfsSUMMary==cur_cf & snrSUMMary==0 & splSUMMARY==cur_spl)==1);
            %             if ~isempty(SummaryVar(check0dBInd).h_revcor)
            %                 temp(i).h_revcor=SummaryVar(check0dBInd).h_revcor;
            %             else
            %                 ind2use=find((extractfield(spike_data,'CF') == cur_cf & extractfield(spike_data,'SPL') == cur_spl & extractfield(spike_data,'SNR') == 0) ==1 );
            %                 if length(ind2use)>1
            %                    ind2use=ind2use(end);
            %                 end
            %             end
            %         else
            ind2use=i;
            %         end
            h=[];
            %         k2=[];
            nSpikes=0;
            %         nSpikesb=0;
            for polarity_var=1:2
                spike_times=cell2mat(spike_data(ind2use).SpikeTrains{2,polarity_var}');
                [pin, Fs]=audioread(StimsFNames{ind2use}{2,polarity_var}{1});
                pin=pin*10^(-cur_attn/20);
                
                [h1,nSpikes1]=find_revcor(pin,spike_times,Fs,window_dur);
                h=[h,h1];
                nSpikes=nSpikes+nSpikes1;
                
                %             [h2,nSpikes2]=find_wiener_kernel(pin,spike_times,Fs,window_dur,1);
                %             [u,~,~]=svd(h2);
                %             k2fsv=u(:,1);
                %             k2=[k2, k2fsv];
                %             nSpikesb=nSpikesb+nSpikes2;
            end
            h_revcor=flipud(mean(h,2));
            
            [SNRenvACST, BF_revcor]=analyze_revcor(StimsFNames{i},h_revcor,cur_attn);
            
            temp(i).h_revcor=h_revcor;
            temp(i).nSpikes4REVCOR=nSpikes;
            temp(i).SNRenvACST=SNRenvACST;
            temp(i).BF_revcor=BF_revcor;
        else
            temp(i).h_revcor=nan;
            temp(i).nSpikes4REVCOR=nan;
            temp(i).SNRenvACST=nan;
            temp(i).BF_revcor=nan;
        end
        %     fprintf('Summary: %d/%d\n',i,length(spike_data));
    end
end

for i=1:length(spike_data)
    if ~isnan(ind2updateInSummaryVar(i))
        
        %    SummaryVar(ind2updateInSummaryVar(i))=temp(i);
        SummaryVar(ind2updateInSummaryVar(i)).h_revcor=temp(i).h_revcor;
        SummaryVar(ind2updateInSummaryVar(i)).nSpikes4REVCOR=temp(i).nSpikes4REVCOR;
        SummaryVar(ind2updateInSummaryVar(i)).SNRenvACST=temp(i).SNRenvACST;
        SummaryVar(ind2updateInSummaryVar(i)).BF_revcor=temp(i).BF_revcor;
    end
end