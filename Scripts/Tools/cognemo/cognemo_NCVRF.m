function cognemo_NCVRF(X,T,options)
%% Preamble
%{

This function computes the generalization error BAC for RF classifiers with
feature selection by t-test, for the dataset group X with labels T.

INPUTS
    X := vectorized datasets [N_o x N_v]
    T := binary class labels [N_o x 1]
    options:
        options.KO    := number of folds for 'outer' CV structure
        options.KI    := number of folds for 'inner' CV structure
        options.seedO := seed for RNG for 'outer' CV partition
        options.seedI := seed for RNG for 'inner' CV partition
        options.N_f   := array of hyperparameter values: number of feats
        options.R_ttf := array of hyperparameter values: ratio,
        trees-to-feats
OUTPUTS
    Inner:
        I_KEEP: indices of connections after filtering (in numeric form, to save space)
        I_TR: indices of training sets for each inner fold
        I_TE: indices of test sets for each inner fold
        % [DELETED] MDL: RF models for each inner fold
        T_PR: predicted labels for each inner fold
        SC: score for label '0' for each inner fold
        T_TE: actual labels for test set for each inner fold
        PM: performance measure values for each inner fold
        PMav: average performance measure values across the inner folds
        OC: optimization criterion values for each inner fold
        OCav: optimization criterion values across the inner folds
        M_conf: confidence matrices for each inner fold
    Outer:
        I_KEEP: indices of connections after filtering (in numeric form, to save space)
        I_TR: indices of training sets for each outer fold
        I_TE: indices of test sets for each outer fold
        N_f: best 'number of features' hyperparameter value for each outer fold
        R_ttf: best 'ratio of trees-to-features' hyperparameter value for each outer fold
        T_PR: predicted labels for each outer fold
        SC: score for label '0' for each outer fold
        T_TE: actual labels for test set for each outer fold
        PM: performance measure values for each outer fold
        PMav: average performance measure values across the outer folds
        OC: optimization criterion value for each outer fold
        OCav: average optimization criterion value across the outer folds
        M_conf: confidence matrices for each outer fold
    
%}
%% Unpack options

% CV structure options
KO = options.KO; seedO = options.seedO;
KI = options.KI; seedI = options.seedI;
N_v = size(X,2);

% RF model hyperparameter options
N_f   = options.N_f;    % hyperparameter array: number of features to select
P     = length(N_f);
R_ttf = options.R_ttf;  % hyperparameter array: ratio of trees-to-features
Q     = length(R_ttf);

% performance measure options
pm_labels = ["ACC","SE","SP","BAC","PPV","NPV","AUC"];
if isfield(options,'pm_labels')
    pm_labels = options.pm_labels;
end
N_pm      = length(pm_labels);
oc_label  = "BAC";
if isfield(options,'oc_label')
    oc_label  = options.oc_label;           % optimization criterion
end
[~,i_oc]  = ismember(oc_label,pm_labels);   % index of oc in pm labels

% for cognemo_classrates
croptions.pm_labels = pm_labels;
croptions.posclass  = options.posclass;

%% Save destinations

% INNER:
    % indices of features kept after filtering
    out.I.I_KEEP = cell([P,Q,KI,KO]);
    % CV partition info
    out.I.I_TR = cell([P,Q,KI,KO]);
    out.I.I_TE = cell([P,Q,KI,KO]);
    % Model, predictions, scores, true labels
    % out.I.MDL    = cell([P,Q,KI,KO]);
    out.I.T_PR   = cell([P,Q,KI,KO]);
    out.I.SC     = cell([P,Q,KI,KO]);
    out.I.T_TE   = cell([P,Q,KI,KO]);
    % performance measures, optimization criterion, confidence matrix
    out.I.PM     = zeros([P,Q,KI,KO,N_pm]);
    out.I.PMav   = zeros([P,Q,KO,N_pm]);
    out.I.OC     = zeros([P,Q,KI,KO]);
    out.I.OCav   = zeros([P,Q,KO]);
    out.I.M_conf = cell([P,Q,KI,KO]);

% OUTER:
    % indices of features kept after filtering
    out.O.I_KEEP = cell([1,KO]);
    % CV partition info
    out.O.I_TR = cell([1,KO]);
    out.O.I_TE = cell([1,KO]);
    % 'best' hyperparameters
    out.O.N_f    = zeros([1,KO]);
    out.O.R_ttf  = zeros([1,KO]);
    % Model, predictions, scores
    % out.O.MDL    = cell([1,KO]);
    out.O.T_PR   = cell([1,KO]);
    out.O.SC     = cell([1,KO]);
    out.O.T_TE   = cell([1,KO]);
    % performance measures, optimization criterion, confidence matrix
    out.O.PM     = zeros([KO,N_pm]);
    out.O.PMav   = zeros([1,N_pm]);
    out.O.OC     = zeros([1,KO]);
    out.O.OCav   = 0;
    out.O.M_conf = cell([1,KO]);

% directories + filenames for saving outer-fold models
dirname = "NCVRF_" + string(datetime); mkdir(dirname); cd(dirname)
OMDLdirname = "outermdl"; mkdir(OMDLdirname); cd ..
OMDLfname = "outermdl_fold";
fname = char(dirname+".mat");

%% Status message

msg.st    = "| %% Done |";
msg.ko    = "| Outer fold |";
msg.ki    = "| Inner fold |";
msg.N_f   = "|  N_f |";
msg.R_ttf = "| R_ttf |";
msg.time  = "|     time |";
printmsg = msg.st+msg.ko+msg.ki+msg.N_f+msg.R_ttf+msg.time+"\n";
fprintf(printmsg)

%% Outer CV structure

rng(seedO); CO = cvpartition(T,'kFold',KO);

n_step = 0; N_step = KO*KI*P*Q;
tic
% begin outer CV loop
for ko = 1:KO
    %% Status message
    msg.ko = "| " + string(repmat(' ',...
                    [1,(strlength("Outer fold") - ...
                    strlength(string(ko)))])) + ...
                    string(ko) + " |";
    
    %% Partition folds
    O = struct;
    % outer training group
    O.I_tr = logical(training(CO,ko));
    O.T_tr = logical(T(O.I_tr)); O.X_tr = X(O.I_tr,:);
    % outer test group
    O.I_te = logical(test(CO,ko));
    O.T_te = logical(T(O.I_te)); O.X_te = X(O.I_te,:);
    
    %% Inner CV structure
    
    rng(seedI); CI = cvpartition(O.T_tr,'kFold',KI);
    
    % begin inner CV loop
    for ki = 1:KI
        %% Status message
        msg.ki = "| " + string(repmat(' ',...
                        [1,(strlength("Inner fold") - ...
                        strlength(string(ki)))])) + ...
                        string(ki) + " |";
        
        %% Partition folds
        I = struct;
        % inner training group
        I.I_tr = logical(training(CI,ki));
        I.T_tr = logical(O.T_tr(I.I_tr)); I.X_tr = O.X_tr(I.I_tr,:);
        % inner test group
        I.I_te = logical(test(CI,ki));
        I.T_te = logical(O.T_tr(I.I_te)); I.X_te = O.X_tr(I.I_te,:);
        
        % Train RF models with each choice of hyperparameters
        for p = 1:P
            %% Status message
            msg.N_f = "| " + string(repmat(' ',...
                             [1,(strlength(" N_f") - ...
                             strlength(string(N_f(p))))])) + ...
                             string(N_f(p)) + " |";
            
            %% select hyperparameter
            I.P.N_f = N_f(p);
            
            % [FUNCTION] feature selection by t-test
            I.P.X0 = I.X_tr(~logical(I.T_tr),:);
            I.P.X1 = I.X_tr(logical(I.T_tr),:);
            [~,~,~,I.P.tt] = ttest2(I.P.X0,I.P.X1);
            [~,I.P.i_keep] = maxk(abs(I.P.tt.tstat),I.P.N_f);
            
            % dimension reduction
            I.P.X_tr_rd = I.X_tr(:,I.P.i_keep);
            I.P.X_te_rd = I.X_te(:,I.P.i_keep);
                        
            for q = 1:Q
                % Progress and time elapsed at this step
                n_step = n_step + 1; progress = round(100*n_step/N_step,2);
                t2 = toc; time = datestr(t2/(24*60*60),"HH:MM:SS");
                
                %% Status message
                msg.time = "| " + string(repmat(' ',...
                                   [1,(strlength("    time") - ...
                                   strlength(time))])) + ...
                                   string(string(time)) + " |"; clear t2 time
                               
                msg.st = "| " + string(repmat(' ',...
                                [1,(strlength("% Done") - ...
                                strlength(string(progress)))])) + ...
                                string(progress) + " |";
                
                msg.R_ttf = "| " + string(repmat(' ',...
                                   [1,(strlength("R_ttf") - ...
                                   strlength(string(R_ttf(q))))])) + ...
                                   string(R_ttf(q)) + " |"; 
                
                %% select hyperparameter
                I.Q.N_tree = R_ttf(q)*I.P.N_f;
                % fit RF model
                I.Q.mdl = TreeBagger(I.Q.N_tree,I.P.X_tr_rd,I.T_tr,...
                                     'Method','classification',...
                                     'PredictorSelection','curvature',...
                                     'OOBPredictorImportance','on');
                                 
                % acquire predictions from RF model
                [I.Q.T_pr,I.Q.sc] = predict(I.Q.mdl,I.P.X_te_rd);
                I.Q.T_pr = str2double(string(cell2mat(I.Q.T_pr)));
                I.Q.sc = I.Q.sc(:,1); I.Q = rmfield(I.Q,'mdl');
                
                % compute performance measures
                [I.Q.rates,I.Q.m_conf,~] = ...
                    cognemo_classrates(I.Q.sc,I.Q.T_pr,I.T_te,I.I_te,croptions);
                
                %% SAVE INNER
                % indices of features kept after filtering
                out.I.I_KEEP(p,q,ki,ko) = {I.P.i_keep};
                % CV partition info
                out.I.I_TR(p,q,ki,ko) = {I.I_tr};
                out.I.I_TE(p,q,ki,ko) = {I.I_te};
                % Model, predictions, scores, true labels
                % out.I.MDL(p,q,ki,ko)  = {I.Q.mdl};
                out.I.T_PR(p,q,ki,ko) = {I.Q.T_pr};
                out.I.SC(p,q,ki,ko)   = {I.Q.sc};
                out.I.T_TE(p,q,ki,ko) = {I.T_te};
                % performance measures, optimization criterion, confidence matrix
                out.I.PM(p,q,ki,ko,:) = I.Q.rates;
                out.I.OC(p,q,ki,ko)   = I.Q.rates(i_oc);
                out.I.M_conf(p,q,ki,ko) = {I.Q.m_conf};
                
                %% Print status
                printmsg = msg.st+msg.ko+msg.ki+msg.N_f+msg.R_ttf+msg.time+"\n";
                fprintf(printmsg)
                
                %%
                I = rmfield(I,'Q');
            end
            I = rmfield(I,'P');
        end
        clear I
    end
    
    % compute oc averaged over inner folds
    O.OCav = mean(out.I.OC(:,:,:,ko),3);   out.I.OCav(:,:,ko) = O.OCav;
    O.PMav = mean(out.I.PM(:,:,:,ko,:),3); out.I.PMav(:,:,ko,:) = O.PMav;
    
    % find p,q indices of best average BAC
    [~,O.ind] = max(O.OCav(:));
    [O.p,O.q] = ind2sub(size(O.OCav),O.ind);
    
    % select best hyperparameter
    O.N_f    = N_f(O.p);
    O.N_tree = R_ttf(O.q)*O.N_f;
    
    % [FUNCTION] feature selection by t-test with best hyperparameter
    O.X0 = O.X_tr(~logical(O.T_tr),:);
    O.X1 = O.X_tr(logical(O.T_tr),:);
    [~,~,~,O.tt] = ttest2(O.X0,O.X1);
    [~,O.i_keep] = maxk(abs(O.tt.tstat),O.N_f);
    
    % dimension reduction
    O.X_tr_rd = O.X_tr(:,O.i_keep); O.X_te_rd = O.X_te(:,O.i_keep);
    
    % fit RF model
    Omdl = TreeBagger(O.N_tree,O.X_tr_rd,O.T_tr,...
                       'Method','classification',...
                       'PredictorSelection','curvature',...
                       'OOBPredictorImportance','on');    
    % acquire predictions from RF model
    [O.T_pr,O.sc] = predict(Omdl,O.X_te_rd);
    O.T_pr = str2double(string(cell2mat(O.T_pr)));
    O.sc = O.sc(:,1);
    
    % save model to file
    cd(dirname); cd(OMDLdirname)
    O.mdlfname = char(OMDLfname+string(ko)+".mat"); save(O.mdlfname,'Omdl');
    cd ../..; clear Omdl
    
    % compute performance measures
    [O.rates,O.m_conf,~] = ...
        cognemo_classrates(O.sc,O.T_pr,O.T_te,O.I_te,croptions);
    
    %% SAVE OUTER
        % indices of features kept after filtering
        out.O.I_KEEP(ko) = {O.i_keep};
        % CV partition info
        out.O.I_TR(ko) = {O.I_tr};
        out.O.I_TE(ko) = {O.I_te};
        % 'best' hyperparameters
        out.O.N_f(ko)   = O.N_f;
        out.O.R_ttf(ko) = R_ttf(O.q);
        % Model, predictions, scores
        % out.O.MDL(ko)  = {O.mdl};
        out.O.T_PR(ko) = {O.T_pr};
        out.O.SC(ko)   = {O.sc};
        out.O.T_TE(ko) = {O.T_te};
        % performance measures, optimization criterion, confidence matrix
        out.O.PM(ko,:)   = O.rates;
        out.O.OC(ko)     = O.rates(i_oc);
        out.O.M_conf(ko) = {O.m_conf};
        
    clear O
end

out.pm_labels = pm_labels;
out.oc_label = oc_label;
out.O.PMav   = mean(out.O.PM,1);
out.O.OCav   = mean(out.O.OC);

fprintf("The generalization error is "+string(out.O.OCav)+"\n")

%% Save output

cd(dirname); save(fname,'out','options'); cd ..

end
