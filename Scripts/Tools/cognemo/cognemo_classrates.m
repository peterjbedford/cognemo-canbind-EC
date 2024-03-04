function [rates,m_conf,t_cont] = cognemo_classrates(T_pr,T_te,idx_te,options)
%% Preamble
%{
This function takes a predicted set of targets and the true (test) set of
targets and computes several performance measures based on classification
hits and misses.
---------------------------------------------------------------------------
INPUTS
---------------------------------------------------------------------------
T_pr :=     predicted targets
T_te :=     actual targets ("te" ~ "test")
options.pm_labels:= string array containing labels for the 
performance measures; may include any combination of:
- "ACC":=    accuracy
- "SE":=     sensitivity
- "SP":=     specificity
- "BAC":=    balanced accuracy
- "PPV":=    positive predictive value
- "NPV":=    negative predictive value
- "AUC":=    area under curve
---------------------------------------------------------------------------
OUTPUTS
---------------------------------------------------------------------------
rates :=    vector containing performance measures
m_conf :=   2-by-2 confusion matrix

%}
%% Unpack defaults

posclass = 0;
if isfield(options,'posclass')
    posclass = options.posclass;
end

pm_labels = ["ACC","SE","SP","BAC","PPV","NPV","AUC"];
if isfield(options,'pm_labels')
    pm_labels = options.pm_labels;
end

N_pm = length(pm_labels);
n_te = length(idx_te);

rates    = zeros(1,N_pm);
m_conf = zeros(2);
t_cont = zeros(n_te,2);

if posclass == 1
    T_te = ~T_te;
    T_pr = ~T_pr;
end
% P==label for 0; N==label for 1;
P_idx  = find(~T_te);        P  = length(P_idx);   % P
TP_idx = find(~T_pr(P_idx)); TP = length(TP_idx);  % TP
FN_idx = find(T_pr(P_idx));  FN = length(FN_idx);  % FN
N_idx  = find(T_te);         N  = length(N_idx);   % N
TN_idx = find(T_pr(N_idx));  TN = length(TN_idx);  % TN
FP_idx = find(~T_pr(N_idx)); FP = length(FP_idx);  % FP

m_conf(1,1) = TP; % number of class 0 correctly classified
m_conf(1,2) = FN; % number of class 0 incorrectly classified as 1
m_conf(2,1) = FP; % number of class 1 incorrectly classified as 0
m_conf(2,2) = TN; % number of class 1 correctly classified

%{
for i = 1:n_te
    t_cont(i,1) = idx_te(i);
    t_cont(i,2) = (T_te(i)==T_pr(i));
end
%}

n = 1;
% accuracy
if ismember("ACC",pm_labels)
    rates(n) = (TP + TN) / (P + N);
    n = n + 1;
end
% sensitivity
if ismember("SE",pm_labels)
    rates(n) = TP / P;
    n = n + 1;
end
% specificity
if ismember("SP",pm_labels)
    rates(n) = TN / N;
    n = n + 1;
end
% balanced accuracy
if ismember("BAC",pm_labels)
    rates(n) = 0.5*( (TP / P) + (TN / N) );
    n = n + 1;
end
% positive predictive value
if ismember("PPV",pm_labels)
    rates(n) = TP / (TP + FP);
    n = n + 1;
end
% negative predictive value
if ismember("NPV",pm_labels)
    rates(n) = TN / (TN + FN);
    n = n + 1;
end
% area under curve
if ismember("AUC",pm_labels)
    [~,~,~,rates(n)] = perfcurve(T_te,scores,posclass);
    n = n + 1;
end

end