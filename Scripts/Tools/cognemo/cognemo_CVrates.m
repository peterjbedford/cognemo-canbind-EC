function rates = cognemo_CVrates(T_pr,T_TE,IDX_TE,options)
%% Preamble
%{
%}
k = options.k;
N_pm = length(options.pm_labels);
rates = zeros(1,N_pm);
n_te = floor(90/k);
for l = 1:k
    ind_l = (1 + (l-1)*n_te):(l*n_te);
    rates = rates + cognemo_classrates(T_pr(ind_l),...
                                       T_TE(ind_l),...
                                       IDX_TE(ind_l),...
                                       options);
end
rates = rates/k;
end