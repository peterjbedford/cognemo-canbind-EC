function p = cognemo_perm(T_pr,T,IDX_TE,rates,options)
%% Preamble
%{
Performs a form of permutation test. The point is to provide a measure of
how much better than chance the classifier performed. Actual target labels
are randomly reordered, then the predicted labels are compared to give a 
'fake' performance measure (BAC) value, and this value is compared with the
'true' BAC value. The proportion of times the fake BAC is less than the
true BAC represents the likelihood the classifier outperforms chance.
---------------------------------------------------------------------------
INPUTS
---------------------------------------------------------------------------
T_pr:=      predicted target labels for test-group datasets
T:=         actual target labels (not important to be in CV-order; just
            need correct proportion of group 0 labels and group 1 labels)
IDX_TE:=    indices of test-group datasets
rates:=     matrix of 'true' performance measure values (contains BAC)
options:=   settings related to CV partition setup
---------------------------------------------------------------------------
OUTPUTS
---------------------------------------------------------------------------
p:=         scalar 'p' value measuring performance compared with chance
%}
%% unpack values
bac = rates(options.pm_labels=="BAC");
N_perm = options.N_perm;
N_o = length(T_pr);

%% permutation test
permoptions = options; permoptions.pm_labels = ["BAC"];
n_better = 0;
for i = 1:N_perm
    ind_i = randperm(N_o);
    T_i = T(ind_i);
    [bac_i] = cognemo_CVrates(T_pr,...
                              T_i,...
                              IDX_TE,...
                              permoptions);
    n_better = n_better + double(bac_i >= bac);
end
p = (n_better + 1) / (N_perm + 1);
end
