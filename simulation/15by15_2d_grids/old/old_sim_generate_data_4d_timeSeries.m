% sim_generate_data (12/21/2013)
% (12/22/2013) <- major updates! updated tak_prec_connectome2d.m big time,
%                 so this script was modified accordingly
% - generate simulation dataset (consumes huge storage space, but whatever)
% - most of the scripts from t_admm_4d_connectome_simulation3.m
%%
clear all
purge
rand('seed',0)
randn('seed',0)

% dispFig=false;
dispFig=true;

rootdir=fileparts(mfilename('fullpath'));
%% set output options
% name of the simulation setup
setup_name = '_setup_june26_2014';

fsave=false;
% fsave=true;
% return

% testing timeseries not saved for storage reason
outVars1={'mFileName','timeStamp','Xtr','Xts','ytr','yts','sim_opt','DISTR','rcorrmean','COV_AR2D'};
outVars2={'mFileName','timeStamp','Ztr','noise_tr'};

outputPath1 = [rootdir,'/sim_data_dataset',setup_name,'.mat']
outputPath2 = [rootdir,'/sim_data_timeseries',setup_name,'.mat']
%% set control parameters 
%==========================================================================
% Parameters on the matrix normal distribution
%==========================================================================
% number of time points
T = 100; 
% temporal correlation
rtime=0.8;
% spatial correlation
rx = 0.25;
ry = 0.25;
% noise level (std dev of noise)
sig = 1;

ntr = 200; 
nts = 500; 

% 2-d node orientation (X,Y)
X=15;
Y=15;
NSIZE=[X Y];

% create precision matrix (2d) from the (vectorized) matrix normal distribution
DISTR = tak_prec_connectome2d(T, rtime, NSIZE,[rx,ry]);
DISTR.arraysize = [DISTR.AR1D.T, DISTR.AR2D.X, DISTR.AR2D.Y];
DISTR.p = prod(DISTR.arraysize);
% return

%==========================================================================
% Other control parameters for the simulation
%==========================================================================
% number of samples (make it even so there are n/2 (+1 class)
%                                          and n/2 (-1 class) samples)
% ntr = 100; 

% number of testing samples 
% (dont' save the timeseries for this.....this is for storage reason)
% nts = 500; 

%==========================================================================
% Regions of the anomalous nodes that I replace the time series
% - x,y patches mirror each other
%==========================================================================
ix_anom = [3:5];  
iy_anom = [11:13]; 

% ix_anom = [3];  
% iy_anom = [4]; 
%==========================================================================
% pool above into a struct 
%==========================================================================
sim_opt.T=T;
sim_opt.rtime=rtime;
sim_opt.rx=rx;
sim_opt.ry=ry;
sim_opt.sig = sig;
sim_opt.ntr=ntr;
sim_opt.nts=nts;
sim_opt.X=X;
sim_opt.Y=Y;
sim_opt.NSIZE=NSIZE;

sim_opt.ix_anom=ix_anom;
sim_opt.iy_anom=iy_anom;
sim_opt
% return
%%
%==========================================================================
% Assign parameters for simulation
% NSIZE = [nx, ny] (# seeds in x,y direction)
% d = nx*ny (total number of seeds)
% p = # connectomes
%==========================================================================
d=prod(sim_opt.NSIZE);
p=d*(d-1)/2;

sim_opt.d=d;
sim_opt.p=p;

% % index list for the +1 subjects (training & testing)
% idxp_tr = 1:ntr/2;
% idxp_ts = 1:nts/2;
% 
% sim_opt.idxp_tr=idxp_tr;
% sim_opt.idxp_ts=idxp_ts;
% 
% % index list for the -1 subjects (training & testing)
% idxm_tr = ntr/2+1:ntr;
% idxm_ts = nts/2+1:nts;
% 
% sim_opt.idxm_tr=idxm_tr;
% sim_opt.idxm_ts=idxm_ts;

% create labels
ytr = [ones(ntr/2,1); -ones(ntr/2,1)];
yts = [ones(nts/2,1); -ones(nts/2,1)];
%==========================================================================
% Draw training and testing spatio-temporal timeseries samples
%   - Ztr (ntr x d*T), d = # nodes, T = # time points
%   - Zts (nts x d*T) 
% - sampling from the inverse covariance requires to solve x=inv(icov_sqrt')*z...
%   ie, the transpose is needed
% - see sanity_check_prec_ar2d_sample_multipletimes.m to remind myself
%   how sampling from our distribution works 
%   (e.g., why there's a transpose below and all that crap)
%==========================================================================
Ztr=(DISTR.A\randn(DISTR.p, ntr))'; % training data timeseries
Zts=(DISTR.A\randn(DISTR.p, nts))'; % testing  data timeseries

%==========================================================================
% reshape the above samples into (T, X, Y, ntr) shape
% (ie, T index the time point, X index x-coord, Y index y-coord,
%      ntr index the sample count)
% - Note: i round using "reshape" directly on Ztr will fuck things up...
%         so be safe, and use loop =)
%==========================================================================
% reshape each sample to a 2-d image...where each pixel has a time series
Zarraytr=zeros(T,X,Y,ntr);
for i =  1:ntr
    Zarraytr(:,:,:,i)=reshape(Ztr(i,:),DISTR.arraysize);
end
Zarrayts=zeros(T,X,Y,nts);
for i =  1:nts
    Zarrayts(:,:,:,i)=reshape(Zts(i,:),DISTR.arraysize);
end

%%% sanity check timeseries plot: check if the timeseries look "right", ie, 
%%% if it agrees with the spatio-temporal correlation structure you assigned. 
%%% note this is a "clean" signal before noise is added
if dispFig
    tak_plot_neighbortime(Zarraytr(:,:,:,1),8,8)
end

% return
%%%%%%%%%%% sanity check %%%%%%%%%%%%%
% => this tells me i have to use the loop above!!!
% tmp_sig=reshape(Ztr,[DISTR.arraysize,ntr]);
% isequal(Zarraytr, tmp_sig)
% tak_plot_neighbortime(Zarraytr(:,:,:,4),3,5)
% tak_plot_neighbortime(tmp_sig(:,:,:,4),3,5)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Ztr=Zarraytr;
Zts=Zarrayts;
clear Zarraytr Zarrayts

% add noise
noise_tr = sim_opt.sig*randn(size(Ztr));
noise_ts = sim_opt.sig*randn(size(Zts));
Ztr = Ztr + noise_tr;
Zts = Zts + noise_ts;

% this time, we plot the noisy version of the sampled signal
if dispFig
    tak_plot_neighbortime(Ztr(:,:,:,1),8,8)
end
% return
%% create "DS" subjects by manually replacing time series within blocks
%==========================================================================
% create the support image of the defective seeds
%==========================================================================
seed_support=zeros(X,Y);
% replace the "defective seeds" (do this for both the training and testing set)
for i = 1:length(ix_anom)
    for j = 1:length(iy_anom)
        ix = ix_anom(i);
        iy = iy_anom(j);
        seed_support(ix,iy) = 1;
        seed_support(iy,ix) = 1;
        
        % replace the timeseries at the anomalous nodes for the +1 subject group --- Z( T, X, Y, isamp)
        Ztr(:,ix,iy,ytr==+1) = reshape( tak_sample_normali(ntr/2,DISTR.AR1D.ICOV)', [T,1,1,ntr/2] );
        Zts(:,ix,iy,yts==+1) = reshape( tak_sample_normali(nts/2,DISTR.AR1D.ICOV)', [T,1,1,nts/2] );

        % the mirrored coordinates (swap ix,iy)
        Ztr(:,iy,ix,ytr==+1) = reshape( tak_sample_normali(ntr/2,DISTR.AR1D.ICOV)', [T,1,1,ntr/2] );
        Zts(:,iy,ix,yts==+1) = reshape( tak_sample_normali(nts/2,DISTR.AR1D.ICOV)', [T,1,1,nts/2] );
    end
end
% return
%%
%==========================================================================
% create connectome signal from sampled time-series
%==========================================================================
% design matrix (vectorized correlations stacked on the rows)
Xtr=zeros(ntr,p); 
for i=1:ntr % Ztr = (T,X,Y,isamp)
    rcorr_tmp=corr(reshape(Ztr(:,:,:,i),[T,d]));
    Xtr(i,:)=tak_dvec(rcorr_tmp);    
%     %%% sanity checks: again, check if the timeseries "looks correct" %%%    
%     tmp=reshape(Ztr(:,:,:,i),[T,d]);
%     tmp2=reshape(tmp,[T,X,Y]);
%     %(8,8) is the x,y coord...the choice of 8 is to keep consistent with earlier plots
%     tak_plot_neighbortime(tmp2,8,8) 
%     % => Note: looks good! the plot looks identical to the one i plotted above =)
%     keyboard
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

% repeat above for the testing set
Xts=zeros(nts,p); 
for i=1:nts % Zts = (T,X,Y,isamp)
    rcorr_tmp=corr(reshape(Zts(:,:,:,i),[T,d]));
    Xts(i,:)=tak_dvec(rcorr_tmp);
end
%% create support matrix for the anomalous node locations
idx_list_anom_lex = []; % lexicographical index list of the anomalous nodes
for i=1:length(ix_anom)
    for j=1:length(iy_anom)
        ix = ix_anom(i);
        iy = iy_anom(j);
        idx_list_anom_lex = [idx_list_anom_lex, sub2ind(NSIZE,ix,iy)];
        idx_list_anom_lex = [idx_list_anom_lex, sub2ind(NSIZE,iy,ix)];
    end
end
idx_list_anom_lex=sort(idx_list_anom_lex);

suppmat=zeros(d);
suppmat(idx_list_anom_lex,:)=1;
suppmat(:,idx_list_anom_lex)=1;

sim_opt.idx_list_anom=idx_list_anom_lex;
sim_opt.suppmat=suppmat;
% if dispFig
%     figure,imagesc(suppmat),axis('image'),colormap(1-gray)
% end
%% get mean correlations/connectomes
%==========================================================================
% Create mean "noise" correlations & "clean-signal" correlations &
%        mean "real-signal" correlations using testing data
% (testing data used since this will remain at 500 samples)
%==========================================================================
% for training data
Xtr_noise =zeros(ntr,p); 
Xtr_clean =zeros(ntr,p); 
for i=1:ntr
    tmp_corr=corr(reshape(noise_tr(:,:,:,i),[T,d]));
    Xtr_noise(i,:)=tak_dvec(tmp_corr);
    
    tmp_corr=corr(reshape(Ztr(:,:,:,i)-noise_tr(:,:,:,i),[T,d]));
    Xtr_clean(i,:)=tak_dvec(tmp_corr);
end

% repeat for testing
Xts_noise =zeros(nts,p); 
Xts_clean =zeros(nts,p); 
for i=1:nts
    tmp_corr=corr(reshape(noise_ts(:,:,:,i),[T,d]));
    Xts_noise(i,:)=tak_dvec(tmp_corr);
    
    tmp_corr=corr(reshape(Zts(:,:,:,i)-noise_ts(:,:,:,i),[T,d]));
    Xts_clean(i,:)=tak_dvec(tmp_corr);
end

%==========================================================================
% compute mean: training
%==========================================================================
rcorrmean.tr_p = tak_dvecinv(mean(Xtr(ytr==+1,:),1));
rcorrmean.tr_m = tak_dvecinv(mean(Xtr(ytr==-1,:),1));

% mean clean signal (signal before noise got added)
rcorrmean.tr_p_clean = tak_dvecinv(mean(Xtr_clean(ytr==+1,:),1));
rcorrmean.tr_m_clean = tak_dvecinv(mean(Xtr_clean(ytr==-1,:),1));

% mean noise
rcorrmean.tr_p_noise = tak_dvecinv(mean(Xtr_noise(ytr==+1,:),1));
rcorrmean.tr_m_noise = tak_dvecinv(mean(Xtr_noise(ytr==-1,:),1));

if dispFig
    figure,imexp
    subplot(231),imcov(rcorrmean.tr_p)
    subplot(234),imcov(rcorrmean.tr_m)
    subplot(232),imcov(rcorrmean.tr_p_clean)
    subplot(235),imcov(rcorrmean.tr_m_clean)
    subplot(233),imcov(rcorrmean.tr_p_noise)
    subplot(236),imcov(rcorrmean.tr_m_noise)
end
%==========================================================================
% compute mean: testing
%==========================================================================
rcorrmean.ts_p = tak_dvecinv(mean(Xts(yts==+1,:),1));
rcorrmean.ts_m= tak_dvecinv(mean(Xts(yts==-1,:),1));

% mean clean signal (signal before noise got added)
rcorrmean.ts_p_clean = tak_dvecinv(mean(Xts_clean(yts==+1,:),1));
rcorrmean.ts_m_clean= tak_dvecinv(mean(Xts_clean(yts==-1,:),1));

% mean noise
rcorrmean.ts_p_noise = tak_dvecinv(mean(Xts_noise(yts==+1,:),1));
rcorrmean.ts_m_noise= tak_dvecinv(mean(Xts_noise(yts==-1,:),1));

if dispFig
    figure,imexp
    subplot(231),imcov(rcorrmean.ts_p),title('rcorrmean,+1')
    subplot(234),imcov(rcorrmean.ts_m),title('rcorrmean,-1')
    subplot(232),imcov(rcorrmean.ts_p_clean),title('rcorrmeanClean,+1')
    subplot(235),imcov(rcorrmean.ts_m_clean),title('rcorrmeanClean,-1')
    subplot(233),imcov(rcorrmean.ts_p_noise),title('noisemean,+1')
    subplot(236),imcov(rcorrmean.ts_m_noise),title('noisemean,-1')
end
%%
%==========================================================================
% display population covariance of the spatial ar covariance, along with
% the covariance model where i zeroed out the location of the anomalous nodes
%==========================================================================
% return
% purge
COV_AR2D.population = full(inv(DISTR.AR2D.ICOV));
COV_AR2D.defective = COV_AR2D.population;
COV_AR2D.defective(idx_list_anom_lex,:)=0;
COV_AR2D.defective(:,idx_list_anom_lex)=0;

% don't want to zero out the diagonal entries....replace values of 1
for iii=1:length(idx_list_anom_lex)
    idx=idx_list_anom_lex(iii);
    COV_AR2D.defective(idx,idx)=1;
end
if dispFig
    figure,imexp
    subplot(231),imcov(rcorrmean.ts_p_clean),title('rcorrmeanClean, +1')
    subplot(232),imcov(COV_AR2D.defective),title('COV_AR2D.defective')
    
    subplot(234),imcov(rcorrmean.ts_m_clean),title('rcorrmeanClean, -1')
    subplot(235),imcov(COV_AR2D.population),title('COV_AR2D.population')
    
    % empirical difference
    diff_emp = rcorrmean.ts_p -   rcorrmean.ts_m;
    
    % population difference
    diff_pop = COV_AR2D.defective   -   COV_AR2D.population;
    
    subplot(233),imcov(diff_emp)
    subplot(236),imcov(diff_pop)

    figure,imexp
    subplot(131),imcov(diff_emp)
    subplot(132),imcov(diff_pop)
    subplot(133),imcov(abs(diff_emp-diff_pop)), title('|diff_emp-diff_pop|')
end
%% save for figures in my report...
timeStamp=tak_timestamp
mFileName=mfilename
% Whos

if fsave
save(outputPath1,outVars1{:})
save(outputPath2,outVars2{:})
end