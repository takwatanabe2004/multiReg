%% do10cv_admm_gridsearch.m
% (05/11/2014)
%=========================================================================%
% - run admm gridsearch - evaluate accuracy via 10fold-CV
%=========================================================================%
%%
clear
purge

% penalty choice
penalty = 'isoTV'; % {'enet','gnet','flas','isoTV'}
%% setup grid search
% assign grid and the name of the grid-configuration
gridname = 'grid_may11';
lamgrid = 2.^(-16:1:-2);     % L1 regularizer
gamgrid = 2.^(-16:1:-2);  % depends on "penalty"

len_gamgrid=length(gamgrid)
len_lamgrid=length(lamgrid)

mFileName=mfilename;
%% manage output path setting and output variables
outputpath=[get_rootdir,'/admm_scripts/gridsearch/',...
            gridname,'_',penalty,'.mat'];
if exist(outputpath,'file')
    reply = input(['***** File already exists at: ', outputpath, ...
                   '\n***** Overwrite? [y/n] '], 's');
    if ~strcmpi(reply,'y'), error('Program terminated'), end
end       

% list of variables to save
outvars={'lamgrid','gamgrid','accuracy','TPR','TNR','nnz_CVcell',...
    'mean_nnz_array','ypred_list','timeStamp','exit_info','wvec_support', ...
    'idx_gamma','idx_lambda','flagDone',...
    'score_list', 'AUC', 'F1',...
    'mFileName','timeStamp'};
% return
%% load data
load sMRI_design_censor.mat X y

[n,p]=size(X);

%==========================================================================
% load info for the indexing for the 10-fold-CV partitioning
%==========================================================================
load Overall Overall
%% set algorithm options
%==========================================================================
% loss function
%==========================================================================
options.loss='hinge1';
% options.loss='hinge2';
% options.loss='hubhinge';
% options.loss_huber_param=0.2; % <- only needed when using huberized-hinge

%==========================================================================
% augmented lagrangian parameters
%==========================================================================
options.rho=1;

%==========================================================================
% termination criterion
%==========================================================================
options.termin.maxiter = 400;   % <- maximum number of iterations
options.termin.tol = 5e-3;      % <- relative change in the primal variable
options.termin.progress = inf;   % <- display "progress" every k iterations
options.termin.silence = true; % <- display termination condition

%==========================================================================
% information needed for data augmentation and fft tricks
%==========================================================================
if ~strcmpi(penalty,'enet')
    
    % augmentation matrix
    load([get_rootdir,'/data_local/A_matrix_cropped.mat'],'A','nx','ny','nz')
    
    % mask vector
    load([get_rootdir,'/data_local/B_matrix_cropped.mat'],'b')

    options.misc.NSIZE=[nx,ny,nz];
    options.misc.A=A; % <- augmentation matrix
    options.misc.b=b; % <- masking vector

    % preassigning this matrix C is optional: 
    % the admm script will compute it internally if the field doesn't exist
    % (helpful when doing CV, since creating this diffmat has minor overhead)
    C=tak_diffmat(options.misc.NSIZE,1); % <- C'*C has circulant structure!
    options.misc.C=C;
end
%%  initialize stuffs to save for gridsearch
ypred_list= cell(len_gamgrid,len_lamgrid);
accuracy = zeros(len_gamgrid,len_lamgrid);
TPR = zeros(len_gamgrid,len_lamgrid);
TNR = zeros(len_gamgrid,len_lamgrid);

% weight vector support and nnz info
wvec_support = cell(len_gamgrid,len_lamgrid);
nnz_CVcell  = cell(len_gamgrid,len_lamgrid);
mean_nnz_array = zeros(len_gamgrid,len_lamgrid);

% struct containing terminatino info
exit_info = cell(len_gamgrid,len_lamgrid);

score_list = cell(len_gamgrid,len_lamgrid);
AUC  = zeros(len_gamgrid,len_lamgrid);
F1 = zeros(len_gamgrid,len_lamgrid);
% return
%% begin gridsearch!
disp(['Begin gridsearch.  Result will be saved at: ',outputpath])
flagDone=false;

for idx_gamma = 1:len_gamgrid %<- \gamma loop
    fprintf('======== idx_gamma = %2d out of %d (gamma=%4.3e) ========\n',...
            idx_gamma,len_gamgrid,gamgrid(idx_gamma))
    % set \gamma value
    options.gamma =gamgrid(idx_gamma);
    for idx_lambda = 1:len_lamgrid %<- \lambda loop
        fprintf('------ idx_lambda = %2d out of %d (lambda=%4.3e)-----\n',...
                idx_lambda,len_lamgrid,lamgrid(idx_lambda))

        % set \lambda value (L1 penalty weight)
        options.lambda=lamgrid(idx_lambda);   
        
        %==================================================================
        % initialization for the CV
        %==================================================================
        ypredicted = [];
        ytrue      = [];
        nnz_tracker= [];
        score = []; % <- NEW! (02/16/2014)
        wvec_support{idx_gamma,idx_lambda}=false(p,10);
        tic_idxCV=tic;
        %%%%%%%%%%%%%%%%%%%%% BEGIN 10-FOLD-CV %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        for idxCV = 1:10
%             fprintf('***** idxCV = %2d ...%6.3f sec *****\n',idxCV,toc(tic_idxCV))
            %=============================================================%
            % 10-fold-CV data partition
            %=============================================================%
            mask_ts = Overall.CrossValidFold(:,idxCV);
            mask_tr = ~mask_ts;

            Xts = X(mask_ts,:);
            Xtr = X(mask_tr,:);

            yts = y(mask_ts);
            ytr = y(mask_tr);

            %=============================================================%
            % run ADMM
            %=============================================================%
            switch penalty
                case 'enet'
                    output=tak_admm_elasticnet(Xtr,ytr,options);
                case 'gnet'
                    output=tak_admm_graphnet(Xtr,ytr,options);
                case 'flas'
                    output=tak_admm_fusedlasso(Xtr,ytr,options);
                case 'isoTV'
                    output=tak_admm_isotropicTV(Xtr,ytr,options);
            end

            %=============================================================% 
            % prediction on test data
            %=============================================================%
            w=output.v2;
            ypr=SIGN(Xts*w);

            ypredicted = [ypredicted;ypr];
            ytrue      = [ytrue; yts];
            
            % score ... used to compute the auc after the LOO 
            score      = [score; Xts*w];

            %=============================================================%
            % weight vector support and nnz info
            %=============================================================%
            wvec_support{idx_gamma,idx_lambda}(:,idxCV)=(w~=0);
            nnz_tracker = [nnz_tracker; nnz(w)];

            %==============================================================
            % terminatino condition info
            %==============================================================
            exit_info{idx_gamma,idx_lambda}.time(idxCV)     = output.time;
            exit_info{idx_gamma,idx_lambda}.k(idxCV)        = output.k;
            exit_info{idx_gamma,idx_lambda}.reltol(idxCV)   = output.rel_change;
        end % <- idxCV
        fprintf('***** idxCV = %2d ...%6.3f sec *****\n',idxCV,toc(tic_idxCV))
        performance=tak_binary_classification_summary(ypredicted,ytrue);
        ypred_list{idx_gamma,idx_lambda}=ypredicted;
        accuracy(idx_gamma,idx_lambda)=performance.accuracy;
        TPR(idx_gamma,idx_lambda)=performance.TPR;
        TNR(idx_gamma,idx_lambda)=performance.TNR;
        
        nnz_CVcell{idx_gamma,idx_lambda}  = nnz_tracker;
        mean_nnz_array(idx_gamma,idx_lambda) = mean(nnz_tracker);
        
        % prediction "scores" (score=X*w)
        score_list{idx_gamma,idx_lambda} = score;
        
        % get auc value
        [trash1,trash2,trash3,AUCval]=perfcurve(ytrue,score,+1);
        AUC(idx_gamma,idx_lambda)=AUCval;
        
        % get F1-score
        F1(idx_gamma,idx_lambda)=performance.F1;
        
        %==================================================================
        % intermediate progress saving
        %==================================================================
        timeStamp=tak_timestamp;
        save(outputpath,outvars{:})
    end % <- idx_lambda
end %<- idx_gamma
%% save
if isfield(options,'misc');
    % for compcat disc storage space, remove field "misc", which contains
    % the edge info for the spatial penalty
    options=rmfield(options,'misc');
end

% pad on variable "options" onto the output list
outvars = [outvars, 'options'];

flagDone=true;
timeStamp=tak_timestamp
save(outputpath,outvars{:})
%% make a figure popout to notify completion =)
figure,imexp