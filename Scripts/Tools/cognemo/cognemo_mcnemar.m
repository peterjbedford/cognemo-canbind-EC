function cognemo_mcnemar(Tc1,Ti1,Tc2,Ti2)
%% Preamble
%{
%}
%%

% get number of disagree-incorrect 1, disagree-incorrect 2,
% agree incorrect:
Nci = length(intersect(find(Tc1),find(Ti2)));
Nic = length(intersect(find(Ti1),find(Tc2)));
Ncc = length(intersect(find(Tc1),find(Tc2)));
Nii = length(intersect(find(Ti1),find(Ti2)));

T_cont = [ Ncc, Nci ;
           Nic, Nii ];

% get test statistic
mcnemar(T_cont)
end