function E_te = cognemo_cctest(X_te,V_te,beta,type)
%% Preamble
%{
This function performs the multivariate regression of the covariates V from
the connectivity data X.
---------------------------------------------------------------------------
INPUTS
---------------------------------------------------------------------------
X_te:=      N_o-by-N_c matrix where N_o is the number of observations and
            N_c is the number of variables. In this case the variables are
            the connection strength parameters.
V_te:=      N_o-by-N_v matrix where N_o is the number of observations and
            N_v is the number of variables to be regressed out.
beta:=      beta parameters from regression on training data.
type:=      "const" adds a row of ones to V so that a constant term (or
            'intercept') can be removed from the data. Otherwise, only the
            betas for the variables in V will be calculated.
---------------------------------------------------------------------------
OUTPUTS
---------------------------------------------------------------------------
E_te:=      error leftover from regression.
%}
% add constant term to 'test' portion of V, if indicated
if type == "const"
    [No_te,~] = size(X_te);
    V_te = [ones(No_te,1), V_te];
end
% remove covariate terms from 'test' data
E_te = X_te - V_te*beta;
end