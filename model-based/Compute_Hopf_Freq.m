clear;

NPARCELLS=1000;
Tmax=440;
NCOND=2;
addpath('..\data');

load ts_meditation1000.mat; % Load your data
%%%
% Parameters of the data
TR=2;  % Repetition Time (seconds)

% Bandpass filter settings
fnq=1/(2*TR);                 % Nyquist frequency
flp = 0.008;                    % lowpass frequency of filter (Hz)
fhi = 0.08;                    % highpass
Wn=[flp/fnq fhi/fnq];         % butterworth bandpass non-dimensional frequency
k=2;                          % 2nd order butterworth filter
[bfilt,afilt]=butter(k,Wn);   % construct the filter
Isubdiag = find(tril(ones(NPARCELLS),-1));

%%%%


for cond=1:NCOND
    xs=tseries(:,cond);
    NSUB=size(find(~cellfun(@isempty,xs)),1);
    fce1=zeros(NSUB,NPARCELLS,NPARCELLS);
    tss=zeros(NPARCELLS,Tmax);
    TT=Tmax;
    Ts = TT*TR;
    freq = (0:TT/2-1)/Ts;
    nfreqs=length(freq);
    PowSpect=zeros(nfreqs,NPARCELLS,NSUB);
    for sub=1:NSUB
        sub
        ts=xs{sub,1};
        [Ns, Tmax]=size(ts);
        TT=Tmax;
        Ts = TT*TR;
        freq = (0:TT/2-1)/Ts;
        nfreqs=length(freq);
        
        for seed=1:NPARCELLS
            x=detrend(ts(seed,:)-mean(ts(seed,:)));
            tss(seed,:)=filtfilt(bfilt,afilt,x);
            pw = abs(fft(tss(seed,:)));
            PowSpect(:,seed,sub) = pw(1:floor(TT/2)).^2/(TT/TR);
        end
        fce1(sub,:,:)=corrcoef(tss(1:NPARCELLS,:)','rows','pairwise');
    end
    
    fce=squeeze(mean(fce1));
    
    Power_Areas=squeeze(mean(PowSpect,3));
    for seed=1:NPARCELLS
        Power_Areas(:,seed)=gaussfilt(freq,Power_Areas(:,seed)',0.01);
    end
    
    [maxpowdata,index]=max(Power_Areas);
    f_diff = freq(index);
    f_diff(find(f_diff==0))=mean(f_diff(find(f_diff~=0)));
    
    save (sprintf('results_f_diff_fce_cond%d.mat', cond),...
        'f_diff', 'fce');
end

