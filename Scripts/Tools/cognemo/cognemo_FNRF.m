function cognemo_FNRF(X,T,options)
%% Preamble
%{  
%}
%% Unpack options

% Important numbers
[N_o,N_c] = size(X); N_r = sqrt(N_c);

% CV structure options
K = options.K; seed = options.seed;

% RF model hyperparameter options
mdlseed = 1;
if isfield(options,'mdlseed')
    mdlseed = options.mdlseed;
end
R_ttf = options.R_ttf;  % ratio of trees-to-features

% performance measure options
pm.labels = ["ACC","SE","SP","BAC","PPV","NPV","AUC"];
if isfield(options,'pm_labels')
    pm.labels = options.pm_labels;
end
pm.N = length(pm.labels);
oc.label  = "BAC";
if isfield(options,'oc.label')
    oc.label  = options.oc_label;           % optimization criterion
end
[~,oc.I]  = ismember(oc.label,pm.labels);   % index of oc in pm labels

% feature importance
imp = 'off';
if isfield(options,'imp')
    imp = options.imp;
end

% covariate correction
cc = 0;
if isfield(options,'cc')
    V = options.cc.V;
end

% for cognemo_classrates
croptions.pm_labels = pm.labels;
croptions.posclass  = options.posclass;

%% Determine networks

fnlabel = options.fnlabel;  % functional network labels for each region

if isfield(options,'fnincl')
    fnincl  = options.fnincl;   % list of networks to compare
else
    fnincl = unique(fnlabel,'stable');
end
fs.N = length(fnincl);
fs.I = zeros(fs.N,N_c); % will contain indices for connections in functional networks
for n = 1:fs.N
    fs.I_r = find(fnlabel==fnincl{n});
    fs.I_c = zeros(N_r); fs.I_c(fs.I_r,fs.I_r) = 1;
    fs.I_c = reshape(fs.I_c,[1,N_c]);
    fs.I(n,:) = fs.I_c; fs = rmfield(fs,{'I_r','I_c'});
end; clear n
out.fs.I = fs.I;

%% Initialize output destinations

% for CV partition info
out.I_TR = zeros(N_o,K);
out.I_TE = zeros(N_o,K);
out.T_TE = zeros(N_o,K);

% for predictions, scores and performance
out.T_PR = zeros(N_o,K,fs.N);
out.PP   = zeros(N_o,K,fs.N);
out.pm.PM     = zeros([pm.N,K,fs.N]);
out.oc.OC     = zeros([K,fs.N]);
out.pm.M_conf = cell([1,K,fs.N]);

% for feature importance
out.imp.IMP = zeros(N_c,K,fs.N);

%% directories + filenames for saving outer-fold models

dirname = "FNRF_" + string(datetime);
cd('Output'); mkdir(dirname); cd(dirname)
MDLdirname = "mdl"; mkdir(MDLdirname);
MDLfname = "mdl_fold";
fname = char(dirname+".mat");

%% Status message

msg.st    = "|  %% Done |";
msg.k     = "|     Fold |";
msg.time  = "|     time |";
msg.fs_N  = "| Network                        |";
msg.time  = "|     time |";
printmsg = msg.st+msg.k+msg.fs_N+msg.time+"\n";
printmsg_orig = printmsg;
fprintf(printmsg)

%% Outer CV structure

rng(seed); C = cvpartition(T,'kFold',K);

n_step = 0; N_step = K*fs.N;
tic

if ~isfield(options,'cc')
    %% NO CC: begin outer CV loop
    for k = 1:K
        %% Status message
        msg.k = "| " + string(repmat(' ',...
                        [1,(strlength("    Fold") - ...
                        strlength(string(k)))])) + ...
                        string(k) + " |";

        %% Partition folds
        O = struct;
        % training group
        O.tr.I = logical(training(C,k));
        O.tr.T = logical(T(O.tr.I)); O.tr.X = X(O.tr.I,:);
        % test group
        O.te.I = logical(test(C,k));
        O.te.T = logical(T(O.te.I)); O.te.X = X(O.te.I,:);

        % Iterate over feature sets
        for n = 1:fs.N

            % dimension reduction
            O.fs.tr.X_rd = O.tr.X(:,logical(fs.I(n,:)));
            O.fs.te.X_rd = O.te.X(:,logical(fs.I(n,:)));

            % number of trees
            O.fs.N_tree = R_ttf*length(find(fs.I(n,:)));

            % fit RF model
            rng(mdlseed);
            mdl = TreeBagger(O.fs.N_tree,O.fs.tr.X_rd,O.tr.T,...
                                 'Method','classification',...
                                 'PredictorSelection','curvature',...
                                 'OOBPredictorImportance',imp);

            % save model to file
            cd(MDLdirname)
            O.fs.mdlfname = char(MDLfname+string(k)+"_"+string(n)+".mat"); save(O.fs.mdlfname,'mdl');
            cd ..;

            % acquire predictions from RF model
            [O.fs.T_pr,O.fs.pp] = predict(mdl,O.fs.te.X_rd);
            O.fs.imp = mdl.OOBPermutedPredictorDeltaError; clear mdl
            O.fs.T_pr = str2double(string(cell2mat(O.fs.T_pr)));
            O.fs.pp = O.fs.pp(:,1);

            % compute performance measures
            [O.fs.rates,O.fs.m_conf,~] = ...
                cognemo_classrates(O.fs.pp,O.fs.T_pr,O.te.T,O.te.I,croptions);

            out.T_PR(O.te.I,k,n) = O.fs.T_pr;
            out.PP(O.te.I,k,n)   = O.fs.pp;
            out.pm.PM(:,k,n)     = O.fs.rates;
            out.oc.OC(k,n)       = O.fs.rates(oc.I);
            out.M_conf(1,k,n)    = {O.fs.m_conf};

            % feature importance
            out.imp.IMP(logical(fs.I(n,:)),k,n) = O.fs.imp;

            %% Status message

            % Progress and time elapsed at this step
            n_step = n_step + 1; progress = round(100*n_step/N_step,2);
            t2 = toc; time = datestr(t2/(24*60*60),"HH:MM:SS");

            % Construct new message
            msg.time = "| " + string(repmat(' ',...
                               [1,(strlength("    time") - ...
                               strlength(time))])) + ...
                               string(string(time)) + " |"; clear t2 time

            msg.st = "| " + string(repmat(' ',...
                            [1,(strlength(" % Done") - ...
                            strlength(string(progress)))])) + ...
                            string(progress) + " |";

            msg.fs_N = "| " + fnincl{n} + ...
                              string(repmat(' ',...
                              [1,(strlength("Network                       ") - ...
                              strlength(fnincl{n}))])) + " |";

            % Print status message
            clc
            fprintf(printmsg_orig)
            printmsg = msg.st+msg.k+msg.fs_N+msg.time+"\n";
            fprintf(printmsg)

            %%
            O = rmfield(O,'fs');

        end

        % Save CV partition info
        out.I_TR(O.tr.I,k) = 1;
        out.I_TE(O.te.I,k) = 1;
        out.T_TE(O.te.I,k) = O.te.T;

        clear O

    end
else
    %% COVARIATE CORRECTION
    for k = 1:K
        %% Status message
        msg.k = "| " + string(repmat(' ',...
                        [1,(strlength("    Fold") - ...
                        strlength(string(k)))])) + ...
                        string(k) + " |";

        %% Partition folds
        O = struct;
        % training group
        O.tr.I = logical(training(C,k));
        O.tr.T = logical(T(O.tr.I)); O.tr.X = X(O.tr.I,:);
        % test group
        O.te.I = logical(test(C,k));
        O.te.T = logical(T(O.te.I)); O.te.X = X(O.te.I,:);
        % for covariate correction
        O.tr.V = V(O.tr.I,:); O.te.V = V(O.te.I,:);

        % Iterate over feature sets
        for n = 1:fs.N

            % dimension reduction
            O.fs.tr.X_rd = O.tr.X(:,logical(fs.I(n,:)));
            O.fs.te.X_rd = O.te.X(:,logical(fs.I(n,:)));

            % number of trees
            O.fs.N_tree = R_ttf*length(find(fs.I(n,:)));

            % obtain residuals
            [beta,O.fs.tr.E_rd] = ...
                            cognemo_cctrain(O.fs.tr.X_rd,O.tr.V,'noconst');
            O.fs.te.E_rd        = ...
                            cognemo_cctest(O.fs.te.X_rd,O.te.V,beta,'noconst');

            % fit RF model
            rng(mdlseed);
            mdl = TreeBagger(O.fs.N_tree,O.fs.tr.E_rd,O.tr.T,...
                                 'Method','classification',...
                                 'PredictorSelection','curvature',...
                                 'OOBPredictorImportance',imp);

            % save model to file
            cd(MDLdirname)
            O.fs.mdlfname = char(MDLfname+string(k)+"_"+string(n)+".mat"); save(O.fs.mdlfname,'mdl');
            cd ..;

            % acquire predictions from RF model
            [O.fs.T_pr,O.fs.pp] = predict(mdl,O.fs.te.E_rd);
            O.fs.imp = mdl.OOBPermutedPredictorDeltaError; clear mdl
            O.fs.T_pr = str2double(string(cell2mat(O.fs.T_pr)));
            O.fs.pp = O.fs.pp(:,1);

            % compute performance measures
            [O.fs.rates,O.fs.m_conf,~] = ...
                cognemo_classrates(O.fs.pp,O.fs.T_pr,O.te.T,O.te.I,croptions);

            out.T_PR(O.te.I,k,n) = O.fs.T_pr;
            out.PP(O.te.I,k,n)   = O.fs.pp;
            out.pm.PM(:,k,n)     = O.fs.rates;
            out.oc.OC(k,n)       = O.fs.rates(oc.I);
            out.M_conf(1,k,n)    = {O.fs.m_conf};

            % feature importance
            out.imp.IMP(logical(fs.I(n,:)),k,n) = O.fs.imp;

            %% Status message

            % Progress and time elapsed at this step
            n_step = n_step + 1; progress = round(100*n_step/N_step,2);
            t2 = toc; time = datestr(t2/(24*60*60),"HH:MM:SS");

            % Construct new message
            msg.time = "| " + string(repmat(' ',...
                               [1,(strlength("    time") - ...
                               strlength(time))])) + ...
                               string(string(time)) + " |"; clear t2 time

            msg.st = "| " + string(repmat(' ',...
                            [1,(strlength(" % Done") - ...
                            strlength(string(progress)))])) + ...
                            string(progress) + " |";

            msg.fs_N = "| " + fnincl{n} + ...
                              string(repmat(' ',...
                              [1,(strlength("Network                       ") - ...
                              strlength(fnincl{n}))])) + " |";

            % Print status message
            clc
            fprintf(printmsg_orig)
            printmsg = msg.st+msg.k+msg.fs_N+msg.time+"\n";
            fprintf(printmsg)

            %%
            O = rmfield(O,'fs');

        end

        % Save CV partition info
        out.I_TR(O.tr.I,k) = 1;
        out.I_TE(O.te.I,k) = 1;
        out.T_TE(O.te.I,k) = O.te.T;

        clear O

    end

out.pm.labels = pm.labels;
out.oc.label  = oc.label;
out.pm.mean   = mean(out.pm.PM,2); out.pm.std = std(out.pm.PM,[],2);
out.oc.mean   = mean(out.oc.OC,1); out.pm.std = std(out.oc.OC,[],1);
out.fnincl    = fnincl;
out.imp.mean  = mean(out.imp.IMP,2); out.imp.std = std(out.imp.IMP,[],2);

[out.best_OC,indbest] = max(out.oc.mean);
out.best_fs = fnincl(indbest);

%% Save output

save(fname,'out','options'); cd ../..

%% Print main results

fprintf("The best performing network was "+string(out.best_fs)+"\n")
fprintf("The value of the optimization criterion was "+string(out.best_OC)+"\n")

end
