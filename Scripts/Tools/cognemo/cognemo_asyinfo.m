function cognemo_asyinfo(C)
%% Preamble
%{
In each condition:
number of significantly different connections across triu/tril
percent of ...
%}
%%

string0_N = ['N_sigdiff = ' char(string(C.asy.tstat0.sig.N_v_sig)) '\n'];
string0_p = ['percent_sigdiff = ' ...
             char(string(C.asy.tstat0.sig.percent_sig)) '\n'];

string1_N = ['N_sigdiff = ' char(string(C.asy.tstat1.sig.N_v_sig)) '\n'];
string1_p = ['percent_sigdiff = ' ...
             char(string(C.asy.tstat1.sig.percent_sig)) '\n'];
         
stringD_N = ['N_sigdiff = ' char(string(C.asy.tstatD.sig.N_v_sig)) '\n'];
stringD_p = ['percent_sigdiff = ' ...
             char(string(C.asy.tstatD.sig.percent_sig)) '\n'];

%% Print Info

outputstr = [ 'Asymmetry Info about ' C.inputname '\n' ...
              'Condition 0:\n' string0_N string0_p ...
              'Condition 1:\n' string1_N string1_p ...
              'Interaction:\n' stringD_N stringD_p ...
            ];

display(sprintf(outputstr))

end