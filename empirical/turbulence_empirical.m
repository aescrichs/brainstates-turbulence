clear;

addpath('..\data');
addpath('..\utils');

load schaefer_MK.mat;
load SClongrange.mat;
load ts_meditation1000.mat;                 % load your data: cell with NSUBxCONDITIONS,
                                            % in each cell a matrix with NPARCELLSxTmax
% Parameters of the data
NPARCELLS=1000;                             % total regions
Tmax=440;                                   % total timepoints
cond=2;                                     % 2 conditions in my timeseries
xs=tseries(:,cond);                         % cell four columns, one for condition
NSUB=size(find(~cellfun(@isempty,xs)),1);   % total subjects in each condition
TR=2;                                       % Repetition Time (seconds)

NR=400;
NRini=20;
NRfin=80;
LAMBDA=[0.30 0.27 0.24 0.21 0.18 0.15 0.12 0.09 0.06 0.03 0.01];
NLAMBDA=length(LAMBDA);

% Bandpass filter settings
fnq=1/(2*TR);                 % Nyquist frequency
flp = 0.008;                  % lowpass frequency of filter (Hz)
fhi = 0.08;                   % highpass
Wn=[flp/fnq fhi/fnq];         % butterworth bandpass non-dimensional frequency
k=2;                          % 2nd order butterworth filter
[bfilt,afilt]=butter(k,Wn);   % construct the filter


for i=1:NPARCELLS
    for j=1:NPARCELLS
        rr(i,j)=norm(SchaeferCOG(i,:)-SchaeferCOG(j,:));
    end
end
range=max(max(rr));
delta=range/NR;

for i=1:NR
    xrange(i)=delta/2+delta*(i-1);
end

enstrophy=zeros(NLAMBDA,NPARCELLS,Tmax);
enstrophy_su=zeros(NLAMBDA,NPARCELLS,Tmax);
signal_filt=zeros(NPARCELLS,Tmax);
Phases=zeros(NPARCELLS,Tmax);
Phases_su=zeros(NPARCELLS,Tmax);
TransferLambda_sub=zeros(NLAMBDA,NSUB);
Turbulence_sub=zeros(NLAMBDA,NSUB);
InformationCascade_sub=zeros(1,NSUB);
Transfer_sub=zeros(1,NSUB);
TransferLambda_su_sub=zeros(NLAMBDA,NSUB);
Turbulence_su_sub=zeros(NLAMBDA,NSUB);
fcr=zeros(NLAMBDA,NSUB);
fcr_su=zeros(NLAMBDA,NSUB);
fclam=zeros(NLAMBDA,NPARCELLS,NPARCELLS);
fclam_su=zeros(NLAMBDA,NPARCELLS,NPARCELLS);


C1=zeros(NLAMBDA,NPARCELLS,NPARCELLS);
[aux indsca]=min(abs(LAMBDA-lambda));
ilam=1;
for lambda=LAMBDA
    for i=1:NPARCELLS
        for j=1:NPARCELLS
            C1(ilam,i,j)=exp(-lambda*rr(i,j));
        end
    end
    ilam=ilam+1;
end

for sub=1:NSUB
    sub
    ts=tseries{sub,cond};
    
    for seed=1:NPARCELLS
        ts(seed,:)=detrend(ts(seed,:)-mean(ts(seed,:)));
        signal_filt(seed,:) =filtfilt(bfilt,afilt,ts(seed,:));
        Xanalytic = hilbert(demean(signal_filt(seed,:)));
        Phases(seed,:) = angle(Xanalytic);
        Phases_su(seed,:)=Phases(seed,randperm(Tmax));
    end
    
    for i=1:NPARCELLS
        for ilam=1:NLAMBDA
            C1lam=squeeze(C1(ilam,:,:));
            sumphases=nansum(repmat(C1lam(i,:)',1,Tmax).*complex(cos(Phases),sin(Phases)))/nansum(C1lam(i,:));
            enstrophy(ilam,i,:)=abs(sumphases);
            % surrogates
            sumphases=nansum(repmat(C1lam(i,:)',1,Tmax).*complex(cos(Phases_su),sin(Phases_su)))/nansum(C1lam(i,:));
            enstrophy_su(ilam,i,:)=abs(sumphases);
        end
    end
    
    for ilam=1:NLAMBDA
        Turbulence_sub(ilam,sub)=nanstd(squeeze(enstrophy(ilam,:)));    
        Turbulence_su_sub(ilam,sub)=nanstd(enstrophy_su(ilam,:));
    end
    
    for ilam=1:NLAMBDA-1
        [cc pp]=corr(squeeze(enstrophy(ilam+1,:,2:end))',squeeze(enstrophy(ilam,:,1:end-1))');
        TransferLambda_sub(ilam+1,sub)=nanmean(abs(cc(find(pp(:)<0.05))));
        %surrogates
        [cc pp]=corr(squeeze(enstrophy_su(ilam+1,:,2:end))',squeeze(enstrophy_su(ilam,:,1:end-1))');
        TransferLambda_su_sub(ilam+1,sub)=nanmean(abs(cc(find(pp(:)<0.05))));        
    end
    
    InformationCascade_sub(sub)=nanmean(TransferLambda_sub(2:NLAMBDA,sub),1);

    %%% Transfer across space
    
    for ilam=1:NLAMBDA
        fclam(ilam,:,:)=corrcoef(squeeze(enstrophy(ilam,:,:))');
        fclam_su(ilam,:,:)=corrcoef(squeeze(enstrophy_su(ilam,:)'));
    end
    
    
    for lam=1:NLAMBDA
        numind=zeros(1,NR);
        fcra=zeros(1,NR);
        for i=1:NPARCELLS
            for j=1:NPARCELLS
                r=rr(i,j);
                index=floor(r/delta)+1;
                if index==NR+1
                    index=NR;
                end
                mcc=fclam(lam,i,j);
                if ~isnan(mcc)
                    fcra(index)=fcra(index)+mcc;
                    numind(index)=numind(index)+1;
                end
            end
        end
        %%% Powerlaw a
        
        grandcorrfcn=fcra./numind;
        clear xcoor;
        clear ycoor;
        nn=1;
        for k=NRini:NRfin
            if grandcorrfcn(k)>0
                xcoor(nn)=log(xrange(k));
                ycoor(nn)=log(grandcorrfcn(k)/grandcorrfcn(NRini));
                nn=nn+1;
            end
        end
        linfunc = @(A, x)(A(1)*x+A(2));
        options=optimset('MaxFunEvals',10000,'MaxIter',1000,'Display','off');
        A0=[-1 1];
        [Afit Residual]= lsqcurvefit(linfunc,A0,xcoor,ycoor,[-4 -10],[4 10],options);
        Transfer_sub(lam,sub)=abs(Afit(1));
        
    end
end
InformationCascade_sub=nanmean(TransferLambda_sub,1);
InformationCascade_su_sub=nanmean(TransferLambda_su_sub,1);

TransferLambda_sub(1,:)=NaN;
TransferLambda_su_sub(1,:)=NaN;


save (sprintf('turbu_measurements%d.mat',cond), 'LAMBDA', 'Transfer_sub', 'InformationCascade_sub',... 
    'Turbulence_sub','InformationCascade_su_sub', 'Turbulence_su_sub', 'TransferLambda_su_sub', 'TransferLambda_sub');
