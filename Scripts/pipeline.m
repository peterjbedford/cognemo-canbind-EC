%% CANBIND ANALYSIS--MDD vs HC

% Step 1: Classification

%% Add tools to path

addpath('Tools/circularGraph')
addpath('Tools/cognemo')

%% Import data

cd Data
A_tab = readtable('ATab.csv');
T_tab = readtable('CTL0MDD1.csv');
V_tab = readtable('CANBIND1-Covariates.csv');
cd ..

%% Assign data to arrays

EC.rdata.X = A_tab{:,1:72900};
data.T = logical(T_tab{:,2});

clear A_tab T_tab

%% Covariate data

data.V = zeros([length(data.T),3]);
MADRS = V_tab{:,2}; EDUC = str2double(string(V_tab{:,4}));
SITE = V_tab{:,3}; SITE_int = zeros([length(data.T),1]);
SITE_u = string(unique(SITE));

for i = 1:length(data.T)
    SITE_int(i) = find(SITE{i}==SITE_u);
end

data.V(:,1) = MADRS;
data.V(:,2) = SITE_int;
data.V(:,3) = EDUC;

clear MADRS SITE SITE_u SITE_int EDUC i V_tab

%% Remove index 81---EDUC value = NA

data.T(81) = [];
EC.rdata.X(81,:) = [];
data.V(81,:) = [];

%% Settings & Prep

EC.xset.symmethod = 2;                      % average triu/tril
% t-test
xset.toptions.corr 	= "none"; % apply FDR correction for multiple comparisons
xset.toptions.alpha = 0.05; % alpha for ttests

%% Downsample

DS_seed = 12;
rng(DS_seed)
N_c = length(EC.rdata.X); N_r = sqrt(N_c);

ind0 = find(~data.T); N0 = length(ind0);
ind1 = find(data.T);  N1 = length(ind1);
Nd = min([N0,N1]); N_o = 2*Nd;
data.Td = zeros([N_o,1]);
EC.rdata.Xd = zeros([N_o,N_c]);
if isfield(data,'V')
    data.Vd = zeros([N_o,3]);
end
if N0 > N1
    ind_d = randi(N0,[Nd,1]);
    ind_d = ind0(ind_d);
    data.Td(1:Nd) = data.T(ind_d);
    data.Td(Nd+1:N_o) = data.T(ind1);
    EC.rdata.Xd(1:Nd,:) = EC.rdata.X(ind_d,:);
    EC.rdata.Xd(Nd+1:N_o,:) = EC.rdata.X(ind1,:);
    if isfield(data,'V')
        data.Vd(1:Nd,:) = data.V(ind_d,:);
        data.Vd(Nd+1:N_o,:) = data.V(ind1,:);
    end
else
    ind_d = randi(N1,[Nd,1]);
    ind_d = ind1(ind_d);
    data.Td(1:Nd) = data.T(ind0);
    data.Td(Nd+1:N_o) = data.T(ind_d);
    EC.rdata.Xd(1:Nd,:) = EC.rdata.X(ind0,:);
    EC.rdata.Xd(Nd+1:N_o,:) = EC.rdata.X(ind_d,:);
    if isfield(data,'V')
        data.Vd(1:Nd,:) = data.V(ind0,:);
        data.Vd(Nd+1:N_o,:) = data.V(ind_d,:);
    end
end

clear N0 N1 ind0 ind1 Nd N_o ind_d

%% RF
%{
%% Classification with CV

% EC.class = cognemo_classoutRF(EC,data,xset);

%% Classification w/out CV

N_tree = 10000;
mdl = TreeBagger(N_tree,EC.pdata.Xd,data.Td,...
                'Method','classification',...
                'PredictorSelection','curvature',...
                'OOBPredictorImportance','on');
ooberrormdl = oobError(mdl);
plot(ooberrormdl)

%% Classification w/ feature selection

N_fs = 1000; N_tree_fs = 10000;

imp = mdl.OOBPermutedPredictorDeltaError;
[~,imp_ind] = sort(imp,'descend');

mdl_data = EC.pdata.Xd(:,imp_ind(1:N_fs));
mdl_fs   = TreeBagger(N_tree_fs,mdl_data,data.Td,...
                'Method','classification',...
                'PredictorSelection','curvature',...
                'OOBPrediction','on');

ooberrormdl_fs = oobError(mdl_fs);
plot(ooberrormdl_fs)
%}

%% SVM + PCA in CV
%{
xset.coptions.N_rd = 100;
EC.class = cognemo_classoutSVM(EC,data,xset);

bac = EC.class.ccout.rates(4)
plot(EC.class.ccout.ROCx,EC.class.ccout.ROCy)
%}

%% t-test
%{
EC.ddata.X = EC.pdata.Xd;
EC.ddata.X0 = EC.ddata.X(~logical(data.Td),:);
EC.ddata.X1 = EC.ddata.X(logical(data.Td),:);
EC.rdataorig = EC.rdata; EC.pdataorig = EC.pdata;
EC.rdata = EC.ddata; EC = rmfield(EC,'pdata');
data.Torig = data.T;
data.T = data.Td;
%}
%% 
%{
EC.pdata = cognemo_prepC(EC,data);
EC = cognemo_mass(EC,xset);

%% nested CV

% real options

options.KO    = 5;
options.KI    = 5;
options.seedO = 12;
options.seedI = 34;
options.N_f   = [100, 200, 300, 500, 1000];
options.R_ttf = [3, 5, 10, 20, 50];
options.posclass = 1;

%{
% practice options
options.KO    = 2;
options.KI    = 2;
options.seedO = 12;
options.seedI = 34;
options.N_f   = [50,100];
options.R_ttf = [3,4];
options.posclass = 1;
%}

X = EC.pdata.X; T = data.T;
cognemo_NCVRF(X,T,options);

%% output X and T

save('in.mat','X','T');
%}
%% Network-by-network

%% Load input
% Load data
% load('in.mat')

data.T_orig     = data.T;
data.V_orig     = data.V;
EC.rdata.X_orig = EC.rdata.X;

data.T = data.Td;
data.V = data.Vd;
EC.rdata.X = EC.rdata.Xd;

data = rmfield(data,{'Td','Vd'});
EC.rdata = rmfield(EC.rdata,'Xd');

EC.pdata = cognemo_prepC(EC,data); 

% Import rlabel variables
load('rlabel.mat');
%% Set options
% defaults -> BAC = 0.67 
%{
options.K = 10;
options.seed = 12;
options.R_ttf = 3;
options.posclass = 1;
options.fnlabel = rlabel.fn_unsort;
%}
options.K = 10;
options.R_ttf = 3;
options.posclass = 1;
options.fnlabel = pset.cgset.rlabel.fns.og;
options.imp = 'on';
options.cc.V = data.V;

X = EC.pdata.X; T = data.T;

% vary random seeds

N_seeds = 50;
%% RUN FNRF
%{
for n = 1:50
    options.seed = n;
    cognemo_FNRF(X,T,options)
end
%}
%% Creating grand array of performance measures
%
dmn.i = 6;
dmn.N_c = 3136;
N_perm = 1000;

PM_grand = zeros(N_seeds,7,options.K,15);
FI_grand = zeros(N_seeds,dmn.N_c,options.K);
perm_grand = zeros(N_seeds,15);
ROCx_grand = cell(N_seeds,15);
ROCy_grand = cell(N_seeds,15);
cd Output; cd FNRF_Feb23;

filenames = strsplit(ls);
for i = 1:((length(filenames)-1)/2)
    filename_i = [filenames{2*i-1} ' ' filenames{2*i}];
    cd(filename_i)
    fname = strtrim(ls('*.mat'));
    fname = string(replace(fname,'cd:','/'));
    load(fname)
    
    PM_grand(i,:,:,:) = out.pm.PM;
    dmn_I = find(out.fs.I(dmn.i,:));
    FI_grand(i,:,:) = out.imp.IMP(dmn_I,:,dmn.i);

    % Permutation test--may not use + ROC
    permoptions.pm_labels = out.pm.labels;
    permoptions.k = options.K; permoptions.N_perm = N_perm;
    for j = 1:15
        perm_grand(i,j) = cognemo_perm(out.T_PR(:,:,j),...
                                       data.T,...
                                       out.I_TE,...
                                       out.pm.PM(:,:,j),...
                                       permoptions);
        I_te_ij = find(out.I_TE);
        scores_ij = out.PP(:,:,j); scores_ij = scores_ij(I_te_ij);
        T_te_ij = out.T_TE; T_te_ij = T_te_ij(I_te_ij);
        
        [ROCx_grand{i,j},ROCy_grand{i,j}] = ...
            perfcurve(T_te_ij,scores_ij,0);
    end
    clear(fname)
    cd ..
end
cd ../..
save('FNRF_grand.mat','PM_grand','FI_grand','perm_grand','ROCx_grand','ROCy_grand')
save('data.mat','data')
%}

%% Stacking

%% Load input
%{
% Load data
load('in.mat')
% Import rlabel variables
load('rlabel.mat');
 %% Set options
options.K  = 5;
options.stack = true;
options.KO = 5;
options.KI = 5;
options.seed = 2;
options.seedO = 2;
options.seedI = 3;
options.R_ttf = 3;
options.N_tree = 500;
options.posclass = 1;
options.fn_label = rlabel.fnlabel_unsort;

cognemo_SFNRF(X,T,options)
%}