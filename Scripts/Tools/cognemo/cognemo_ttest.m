function toutput = cognemo_ttest(X0,X1,toptions)
%% Preamble
%{
%}
%%
% Split dataset into observations in condition 1, condition 2
alpha = 0.05;
if isfield(toptions,'alpha')
    alpha = toptions.alpha;
end
[toutput.h,toutput.p,~,stats] = ttest(X0,X1,'alpha',alpha);
toutput.tstat = stats.tstat;
toutput.sd    = stats.sd;
if isfield(toptions,'corr') && (toptions.corr == "FDR")
    % FDR correction for multiple comparisons
    [toutput.hc,toutput.crit_p,~,toutput.adj_p] = fdr_bh(toutput.p);
end
end
