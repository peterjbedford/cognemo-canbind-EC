function T_pr = cognemo_score2class(SCORES,th)
%% Preamble
%{
%}
%% 
ind0 = SCORES>th;
ind1 = SCORES<th;
indt = SCORES==th;

ties = randi([0,1],size(find(indt)));

T_pr = zeros(size(SCORES));
T_pr(ind0) = 0;
T_pr(ind1) = 1;
T_pr(indt) = ties;

end