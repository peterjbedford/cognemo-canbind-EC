function cognemo_tinfo(C)
%% Preamble
%{
number of regions
number of connections
number of features kept for analysis
number of significantly different connections across conditions
percent of ...
%}
%%
string_N_r   = ['Number of regions N_r = ' char(string(sqrt(C.pdata.N_c))) '\n'];
string_N_c   = ['Number of connections N_c = ' char(string(C.pdata.N_c)) '\n'];
string_N_v   = ['Number of features kept for analysis N_v =' ...
                char(string(C.pdata.N_v)) '\n'];
string_N_sig = [...
        'Number of statistically significantly different features N_sig = ' ...
        char(string(C.tstat.sig.N_v_sig)) '\n'];
string_p_sig = [...
       'Percent of statistically significantly different features N_sig = ' ...
       char(string(C.tstat.sig.percent_sig)) '\n'];

%% Print Info

outputstr = [ 'Info about ' C.inputname '\n' ...
              string_N_r string_N_c string_N_v ...
              string_N_sig string_p_sig ...
            ];

display(sprintf(outputstr))

end