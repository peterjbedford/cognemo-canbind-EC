function [beta,E_tr] = cognemo_cctrain(X_tr,V_tr,type)
%% Preamble
%{
This function performs the multivariate regression of the covariates V from
the connectivity data X.
---------------------------------------------------------------------------
INPUTS
---------------------------------------------------------------------------
X:=         N_o-by-N_c matrix where N_o is the number of observations and
            N_c is the number of variables. In this case the variables are
            the connection strength parameters.
V:=         N_o-by-N_v matrix where N_o is the number of observations and
            N_v is the number of variables to be regressed out.
type:=      "const" adds a row of ones to V so that a constant term (or
            'intercept') can be removed from the data. Otherwise, only the
            betas for the variables in V will be calculated.
---------------------------------------------------------------------------
OUTPUTS
---------------------------------------------------------------------------
beta:=      beta parameters from regression.
E:=         error leftover from regression.
%}
%%
if type=="const"
    [No,~] = size(X_tr);
    V_tr = [ones(No, 1), V_tr];
end
    beta = pinv(V_tr)*X_tr;
    E_tr    = X_tr - V_tr*beta;
end
