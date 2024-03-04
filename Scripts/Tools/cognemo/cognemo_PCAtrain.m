function [MX_tr,PX_tr,X_tr_rd] = cognemo_PCAtrain(X_tr,N_rd)
%% Preamble
%{
%}
%%

X_tr = X_tr';
% de-mean
MX_tr = mean(X_tr, 2);
X_tr = X_tr - MX_tr;
% do PCA by svds
[U,~,~] = svds(X_tr,N_rd);

% reduce
PX_tr = U;
X_tr_rd = PX_tr'*X_tr; X_tr_rd = X_tr_rd';

end
