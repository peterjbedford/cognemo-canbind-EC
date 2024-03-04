function cognemo_SFNRF(X,T,options)
%% Preamble
%{
%}
%% Unpack options

% Do stacking?
stack = false;
if isfield(options,'stack')
    stack = options.stack;
end

% Important numbers
[N_o,N_c] = size(X); N_r = sqrt(N_c);

% CV structure options
if ~stack
    KO = options.K; seedO = options.seed;
else
    KO = options.KO; seedO = options.seedO;
    KI = options.KI; seedI = options.seedI;
end

% RF model hyperparameter options
mdlseed = 1;
if isfield(options,'mdlseed')
    mdlseed = options.mdlseed;
end
R_ttf  = options.R_ttf;  % ratio of trees-to-features, for base models
N_tree = options.N_tree;  % number of trees, for meta models

% performance measure options
pm.labels = ["ACC","SE","SP","BAC","PPV","NPV","AUC"];
if isfield(options,'pm_labels')
    pm.labels = options.pm_labels;
end
pm.N      = length(pm.labels);
oc.label  = "BAC";
if isfield(options,'oc.label')
    oc.label  = options.oc_label;           % optimization criterion
end
[~,oc.I]  = ismember(oc.label,pm.labels);   % index of oc in pm labels

% for cognemo_classrates
croptions.pm_labels = pm.labels;
croptions.posclass  = options.posclass;

%% Determine networks

fn.label = options.fn_label;  % functional network labels for each region

if isfield(options,'fn.incl')
    fn.incl = options.fn_incl;   % list of networks to compare
else
    fn.incl = unique(fn.label,'stable');
end
fs.N = length(fn.incl);
fs.I = zeros(fs.N,N_c); % will contain indices for connections in functional networks 
for n = 1:fs.N
    fs.I_r = find(fn.label==fn.incl{n});
    fs.I_c = zeros(N_r); fs.I_c(fs.I_r,fs.I_r) = 1;
    fs.I_c = reshape(fs.I_c,[1,N_c]);
    fs.I(n,:) = fs.I_c; fs = rmfield(fs,{'I_r','I_c'});
end; clear n
out.base.I = fs.I;

%% Initialize output destinations 

% for CV partition info
out.I_TR = zeros(N_o,KO);
out.I_TE = zeros(N_o,KO);
out.T_TE = zeros(N_o,KO);

% for metamodels' predictions, scores and performance
out.meta.T_PR = zeros(N_o,KO);
out.meta.PP   = zeros(N_o,KO);
out.meta.pm.PM  = zeros([pm.N,KO]);
out.meta.oc.OC  = zeros([1,KO]);
out.meta.M_conf = cell([1,KO]);

% for basemodels' predictions, scores and performance
out.base.T_PR = zeros(N_o,KO,fs.N);
out.base.PP   = zeros(N_o,KO,fs.N);
out.base.pm.PM  = zeros([pm.N,KO,fs.N]);
out.base.oc.OC  = zeros([KO,fs.N]);
out.base.M_conf = cell([1,KO,fs.N]);

%% directories + filenames for saving outer-fold models

dirname = "SNVRF_" + string(datetime); mkdir(dirname); cd(dirname)
MMDLdirname = "metamdl"; mkdir(MMDLdirname); cd ..
MMDLfname = "metamdl_fold";
fname = char(dirname+".mat");

%% Status message

msg.st    = "| %% Done |";
msg.ko    = "| Outer fold |";
msg.ki    = "| Inner fold |";
msg.fs_N  = "| Network                        |";
msg.time  = "|     time |";
printmsg = msg.st+msg.ko+msg.ki+msg.fs_N+msg.time+"\n";
printmsg_orig = printmsg;
fprintf(printmsg)

%% Outer CV structure

rng(seedO); CO = cvpartition(T,'kFold',KO);

n_step = 0;
if ~stack
    N_step = KO*fs.N;
else
    N_step = KO * (KI + 1) * fs.N;
end

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
    O.tr.I = logical(training(CO,ko));
    O.tr.T = logical(T(O.tr.I)); O.tr.X = X(O.tr.I,:);
    % outer test group
    O.te.I = logical(test(CO,ko));
    O.te.T = logical(T(O.te.I)); O.te.X = X(O.te.I,:);
    
    if stack
        O.tr.N = length(O.tr.T);
        O.te.N = length(O.te.T);
        % initialize predicted probabilities storage
        O.base.tr.pp = zeros(O.tr.N,fs.N);
        O.base.te.pp = zeros(O.te.N,fs.N);
    end
    
    % iterate over feature sets
    for n = 1:fs.N

        % dimension reduction
        O.fs.tr.X_rd = O.tr.X(:,logical(fs.I(n,:)));
        O.fs.te.X_rd = O.te.X(:,logical(fs.I(n,:)));

        % number of trees
        O.fs.N_tree = R_ttf*length(find(fs.I(n,:)));
        % fit base model
        rng(mdlseed);
        basemdl = TreeBagger(O.fs.N_tree,O.fs.tr.X_rd,O.tr.T,...
                            'Method','classification',...
                            'PredictorSelection','curvature',...
                            'OOBPredictorImportance','off');
                        
        % acquire predicted probabilities for outer test (validation) set
        [O.fs.T_pr,O.fs.pp] = predict(basemdl,O.fs.te.X_rd); clear basemdl
        O.fs.T_pr = str2double(string(cell2mat(O.fs.T_pr)));
        O.fs.pp = O.fs.pp(:,1); 

        % save to storage
        if stack    O.base.te.pp(:,n) = O.fs.pp;    end
        
        % compute performance measures
        [O.fs.rates,O.fs.m_conf,~] = ...
            cognemo_classrates(O.fs.pp,O.fs.T_pr,O.te.T,O.te.I,croptions);

        out.base.T_PR(O.te.I,ko,n) = O.fs.T_pr;
        out.base.PP(O.te.I,ko,n)   = O.fs.pp;
        out.base.pm.PM(:,ko,n)     = O.fs.rates;
        out.base.oc.OC(ko,n)       = O.fs.rates(oc.I);
        out.base.M_conf(1,ko,n)    = {O.fs.m_conf};
        
        %% Status message

        % Progress and time elapsed at this step
        n_step = n_step + 1; progress = round(100*n_step/N_step,2);
        t2 = toc; time = datestr(t2/(24*60*60),"HH:MM:SS");

        % Construct new status message

        msg.time = "| " + string(repmat(' ',...
                           [1,(strlength("    time") - ...
                           strlength(time))])) + ...
                           string(string(time)) + " |"; clear t2 time

        msg.st = "| " + string(repmat(' ',...
                        [1,(strlength("% Done") - ...
                        strlength(string(progress)))])) + ...
                        string(progress) + " |";

        msg.fs_N = "| " + fn.incl{n} + ...
                          string(repmat(' ',...
                          [1,(strlength("Network                       ") - ...
                          strlength(fn.incl{n}))])) + " |";

        % Print status message
        clc
        fprintf(printmsg_orig)
        printmsg = msg.st+msg.ko+msg.ki+msg.fs_N+msg.time+"\n";
        fprintf(printmsg)

        %%
        O = rmfield(O,'fs');

    end
    
    if stack
        %% Inner CV structure

        rng(seedI); CI = cvpartition(O.tr.T,'kFold',KI);

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
            I.tr.I = logical(training(CI,ki));
            I.tr.T = logical(O.tr.T(I.tr.I)); I.tr.X = O.tr.X(I.tr.I,:);
            % inner test group
            I.te.I = logical(test(CI,ki));
            I.te.T = logical(O.tr.T(I.te.I)); I.te.X = O.tr.X(I.te.I,:);
            I.te.N = length(I.te.T);

            % initialize storage
            I.te.pp = zeros(I.te.N,fs.N);

            % begin feature sets loop
            for n = 1:fs.N

                % dimension reduction
                I.fs.tr.X_rd = I.tr.X(:,logical(fs.I(n,:)));
                I.fs.te.X_rd = I.te.X(:,logical(fs.I(n,:)));

                % number of trees
                I.fs.N_tree = R_ttf*length(find(fs.I(n,:)));
                % fit base model
                rng(mdlseed)
                basemdl = TreeBagger(I.fs.N_tree,I.fs.tr.X_rd,I.tr.T,...
                                    'Method','classification',...
                                    'PredictorSelection','curvature',...
                                    'OOBPredictorImportance','off');

                % acquire predicted probabilities for test set
                [~,I.fs.pp] = predict(basemdl,I.fs.te.X_rd); clear basemdl
                I.fs.pp = I.fs.pp(:,1);

                % store predicted probabilities
                I.te.pp(:,n) = I.fs.pp;

                %% Status message

                % Progress and time elapsed at this step
                n_step = n_step + 1; progress = round(100*n_step/N_step,2);
                t2 = toc; time = datestr(t2/(24*60*60),"HH:MM:SS");

                % Construct new status message
                msg.time = "| " + string(repmat(' ',...
                                   [1,(strlength("    time") - ...
                                   strlength(time))])) + ...
                                   string(string(time)) + " |"; clear t2 time

                msg.st = "| " + string(repmat(' ',...
                                [1,(strlength("% Done") - ...
                                strlength(string(progress)))])) + ...
                                string(progress) + " |";

                msg.fs_N = "| " + fn.incl{n} + ...
                                  string(repmat(' ',...
                                  [1,(strlength("Network                       ") - ...
                                  strlength(fn.incl{n}))])) + " |";

                % Print status
                clc
                fprintf(printmsg_orig)
                printmsg = msg.st+msg.ko+msg.ki+msg.fs_N+msg.time+"\n";
                fprintf(printmsg)

                %%
                I = rmfield(I,'fs');

            end

            % save to outer fold storage
            O.base.tr.pp(I.te.I,:) = I.te.pp;

            clear I
        end
        
        % fit meta model
        
        rng(mdlseed);
        metamdl = TreeBagger(N_tree,O.base.tr.pp,O.tr.T,...
                            'Method','classification',...
                            'PredictorSelection','curvature',...
                            'OOBPredictorImportance','off');
        % acquire final predicted probabilities for meta model
        [O.meta.T_pr,O.meta.pp] = predict(metamdl,O.base.te.pp);
        O.meta.T_pr = str2double(string(cell2mat(O.meta.T_pr)));
        O.meta.pp = O.meta.pp(:,1);
        
        %{
        metamdl = fitcdiscr(O.base.tr.pp,O.tr.T);
        % try LDA instead of fitclinear
        
        % acquire final predicted probabilities for meta model
        [O.meta.T_pr,O.meta.pp] = predict(metamdl,O.base.te.pp);
        %O.meta.T_pr = str2double(string(cell2mat(O.meta.T_pr)));
        O.meta.pp = O.meta.pp(:,1);
        %}
        % save model to file
        cd(dirname); cd(MMDLdirname)
        O.mdlfname = char(MMDLfname+string(ko)+".mat"); save(O.mdlfname,'metamdl');
        cd ../..; clear metamdl

        % compute performance measures
        [O.meta.rates,O.meta.m_conf,~] = ...
            cognemo_classrates(O.meta.pp,O.meta.T_pr,O.te.T,O.te.I,croptions);

        % save to output
        out.meta.T_PR(O.te.I,ko,n) = O.meta.T_pr;
        out.meta.PP(O.te.I,ko,n)   = O.meta.pp;
        out.meta.pm.PM(:,ko) = O.meta.rates;
        out.meta.oc.OC(ko)   = O.meta.rates(oc.I);
        out.meta.M_conf(ko)  = {O.meta.m_conf};
    
    end
    
    % Save CV partition info
    out.I_TR(O.tr.I,ko) = 1;
    out.I_TE(O.te.I,ko) = 1;
    out.T_TE(O.te.I,ko) = O.te.T;
    
    clear O 
    
end

out.pm.labels = pm.labels;
out.oc.label  = oc.label;

out.base.pm.mean = mean(out.base.pm.PM,2); out.base.pm.std = std(out.base.pm.PM,[],2);
out.base.oc.mean = mean(out.base.oc.OC,1); out.base.pm.std = std(out.base.oc.OC,[],1);

if stack
    out.meta.pm.mean = mean(out.meta.pm.PM,2); out.meta.pm.std = std(out.meta.pm.PM,[],2);
    out.meta.oc.mean = mean(out.meta.oc.OC);   out.meta.oc.std = std(out.meta.oc.OC);
end

out.fnincl    = fn.incl;

[out.best_OC,indbest] = max(out.base.oc.mean);
out.best_fs = fn.incl(indbest);

%% Save output

cd(dirname); save(fname,'out','options'); cd ..

%% Print main results

fprintf("The best performing network was "+string(out.best_fs)+"\n")
fprintf("The value of the optimization criterion was "+string(out.best_OC)+"\n")

if stack
    fprintf("The generalization error for the metamodels was "+string(out.meta.oc.mean)+"\n")
end

end
