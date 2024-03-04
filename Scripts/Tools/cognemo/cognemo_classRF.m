function [SCORES,IMP,IDX_TE,T_TE] = cognemo_classRF(X,T,options)
%% Preamble
%{
%}
%% unpack options

% Number of CV folds
k = 10;     % DEFAULT
if isfield(options,'k')
    k = options.k;
end

% reduction by t-test?
tt = 0;
if isfield(options,'N_tt')
    tt = 1;
    N_tt = options.N_tt;
end

treestofeats = 1;
if isfield(options,'treestofeats')
    treestofeats = options.treestofeats;
end

%% prep

rng(1234)

% Create CV partition
cvp = cvpartition(T,'kFold',k);

% Create banks for prediction scores and feature importance
[N_o,N_c] = size(X);
SCORES = zeros(1,N_o);
IMP    = zeros(k,N_c);
IDX_TE = zeros(1,N_o);

%% CV loop

idx_l = 0;
for l = 1:k
    fprintf(['Starting fold ' char(string(l)) '\n']);
    
    idx_tr = find(training(cvp,l)); idx_te = find(test(cvp,l));
    X_tr = X(idx_tr,:); X_te = X(idx_te,:);
    T_tr = T(idx_tr);   T_te = T(idx_te);
    
    idx_l = (1:length(idx_te)) + max(idx_l);
    
    % t-test reduction 
    if tt
        ind0 = find(~T_tr); ind1 = find(T_tr);
        No_tt = min([length(ind0),length(ind1)]);
        X0 = X_tr(ind0(1:No_tt),:); X1 = X_tr(ind1(1:No_tt),:);
        [h,~,~,stats] = ttest(X0,X1); N_keep = min(length(find(h)),N_tt);
        [~,ind_keep] = mink(stats.tstat,N_keep);
        X_tr = X_tr(:,ind_keep); X_te = X_te(:,ind_keep);
        options.N_tree = treestofeats*N_keep;
    end

    % CLASSIFICATION
    mdl        = cognemo_RFtrain(X_tr,T_tr,options);
    [~,score0] = cognemo_RFtest(mdl,X_te);

    % EXTRACT IMPORTANT FEATURES
    imp = mdl.OOBPermutedPredictorDeltaError;

    % SAVE
    SCORES(idx_l) = score0;
    if tt  IMP(l,ind_keep) = imp;  else  IMP(l,:) = imp;  end
    IDX_TE(idx_l) = idx_te;
    T_TE(idx_l)   = T_te;

    clear idx_tr idx_te X_tr X_te T_tr T_te ...
          mdl score0 imp
    fprintf(['Done fold ' char(string(l)) '\n'])
end  


end