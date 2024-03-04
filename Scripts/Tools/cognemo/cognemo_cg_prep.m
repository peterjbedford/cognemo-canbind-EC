function Xplot = cognemo_cg_prep(X,rlabel_ind)
%% Preamble
%{
---------------------------------------------------------------------------
INPUTS
---------------------------------------------------------------------------
X must be a vectorized symmetric adjacency matrix
---------------------------------------------------------------------------
OUTPUTS
---------------------------------------------------------------------------
Xplot is a symmetric adjacency matrix in matrix form and reordered to match
plot labels
%}
%%
N_c = length(X); N_r = sqrt(N_c);

X_mat = reshape(X,[N_r,N_r]); 
Xplot = X_mat(rlabel_ind,rlabel_ind); % reorder to match plot labels

end
