function X = cognemo_getFC(data)
%%
%{
This function computes the functional connectivity matrix for a folder
dataloc of datasets, saving the results to separate files in a folder
saveloc.
---------------------------------------------------------------------------
INPUTS
---------------------------------------------------------------------------
data:= data.Yall is a cell of structures, with each struct containing:
            Y.y     := length-Nt BOLD serie,s for Nr regions (NtxNr matrix)
            Y.dt    := scalar value of BOLD time step duration
            Y.name  := optional cell of string labels for regions
            Y.cond  := boolean or 0/1 integer value indicating group 
                       membership or condition
            Y.subj  := subject number
---------------------------------------------------------------------------
OUTPUTS
---------------------------------------------------------------------------
X:=     matrix containing vectorized correlation matrices
%}
%%
Yall = data.Yall;

N_d = size(Yall,2);
N_r = size(Yall{1,1}.y,2);
N_c = N_r*N_r;

X = zeros(N_d,N_c);

for i=1:N_d
    clear yi Si Xi
    
    % Compute FC
    yi = Yall{1,i}.y;
    Si = std(yi);
    Xi = cov(yi)./(Si'*Si);
    
    % Reshape
    X(i,:) = reshape(Xi,[1,N_c]);
end

end