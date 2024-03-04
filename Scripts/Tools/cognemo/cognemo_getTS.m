function [T,S] = cognemo_getTS(Yall)
%% Preamble
%{
Gets condition (T) and subject (S) labels for datasets in Yall.
---------------------------------------------------------------------------
INPUTS
---------------------------------------------------------------------------
Yall:= Yall is a cell of structures, with each struct containing:
            Y.y     := length-Nt BOLD serie,s for Nr regions (NtxNr matrix)
            Y.dt    := scalar value of BOLD time step duration
            Y.name  := optional cell of string labels for regions
            Y.cond  := boolean or 0/1 integer value indicating group 
                       membership or condition
            Y.subj  := subject number
---------------------------------------------------------------------------
OUTPUTS
---------------------------------------------------------------------------
T:=     length-N_o (N_o is number of observations) vector containing
        condition values (collected from the separate Y.cond values)
S:=     length-N_o vector containing subject numbers (indexes of T and S
        match)
%}
%%
N_d = size(Yall,2);

T = zeros(N_d,1); S = zeros(N_d,1);
for i = 1:N_d
    T(i) = Yall{1,i}.cond;
    S(i) = Yall{1,i}.subj;
end

T = logical(T);

end