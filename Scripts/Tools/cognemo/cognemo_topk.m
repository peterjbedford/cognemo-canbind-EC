function [Xtop,Xtop_ind] = cognemo_topk(X,k,method)
%% Preamble
%{
%}
%%
N_c = length(X);
if method == 1
    [~,Xtop_ind] = maxk(abs(X),k);
else
    [~,Xtop_ind] = maxk(X,k);
end
Xtop = zeros(size(X));
Xtop(Xtop_ind) = X(Xtop_ind);
end