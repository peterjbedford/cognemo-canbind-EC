function X_te_rd = cognemo_PCAtest(X_te,MX_tr,PX_tr)
%% Preamble
%{
%}
%%

X_te = X_te';

% de-mean
X_te = X_te - MX_tr;
% reduce
X_te_rd = PX_tr'*X_te;
X_te_rd = X_te_rd';

end
