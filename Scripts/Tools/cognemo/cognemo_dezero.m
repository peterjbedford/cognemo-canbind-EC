function [X_nz,nz_ind,nz_f_ind,z_ind,z_f_ind] = cognemo_dezero(X)
%% Preamble
%{
This function reduces the data matrix X such that any connections that are
zero for all observations are removed. This can be used to reduce the
matrix size before running classification.
---------------------------------------------------------------------------
INPUTS
---------------------------------------------------------------------------
X:=         N_o-by-N_c matrix where N_o is the number of observations and
            N_c is the number of variables. In our case the variables are
            connection strength parameters.
---------------------------------------------------------------------------
OUTPUTS
---------------------------------------------------------------------------
X_dz:=      N_o-by-N_f matrix where N_f is the number of remaining
            variables after removing the all-zero connections.
f_ind:=     Indices of remaining connections.
z_ind:=     Indices of all-zero connections.

%}
%%
nz_ind = find(~all(X==0,1));  % features to keep
nz_f_ind = zeros([1,size(X,2)]); nz_f_ind(nz_ind) = 1;
nz_f_ind = logical(nz_f_ind);
z_ind = find(all(X==0)); % all zero columns
z_f_ind = zeros([1,size(X,2)]); z_f_ind(z_ind) = 1;
z_f_ind = logical(z_f_ind);
X_nz = X(:,nz_ind);

end