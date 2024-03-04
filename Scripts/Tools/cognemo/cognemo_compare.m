function [MDXkeep,SDXkeep,tstatkeep,keep_ind] = cognemo_compare(X0,X1,toptions)
%% Preamble
%{
%}
%% Default n_keep 

n = size(X0,2);
if isfield(toptions,'n')
    n = toptions.n;
end

%% Mean, std

DX = X0 - X1; % (+) implies increased EC in Condition 0 than in Condition 1
MDX = mean(DX,1); SDX = std(DX,1);

%% Significance

toutput = cognemo_ttest(X0,X1,toptions);
% t-values
tstat   = toutput.tstat;
% h=1 if null rejected, h=0 if not
h = toutput.h;
if isfield(toutput,'hc')
    % change to corrected version if indicated
    h = toutput.hc;
end
% choose output to keep
if length(find(h)) <= n
    % choose top values if provided n_keep is less than number significant
    keep_ind = find(h);
else
    % otherwise, must take those that are 'most' significant
    tstath = tstat(h);
    [~,t_sort_ind] = sort(abs(tstath),'descend');
    keep_ind_unsort     = t_sort_ind(1:n);
    [keep_ind,~]        = sort(keep_ind_unsort,'ascend');
end

tstatkeep = tstat(keep_ind);
MDXkeep = MDX(keep_ind); SDXkeep = SDX(keep_ind);

end