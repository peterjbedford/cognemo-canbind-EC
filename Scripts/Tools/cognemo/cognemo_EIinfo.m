function cognemo_EIinfo(C,result)
%% Preamble
%{
%}
%%
if result
    N_c = C.pdata.N_c; N_r = sqrt(N_c);

    % self-connections
    di_ind = logical(reshape(eye(N_r),[1,N_c]));
    X0.Xdi = C.pdata.X0_f(:,di_ind); X0.MXdi = mean(X0.Xdi,'all');
    X1.Xdi = C.pdata.X1_f(:,di_ind); X1.MXdi = mean(X1.Xdi,'all');

    di_string1 = ['a_mean = ' char(string(X0.MXdi)) '\n'];
    di_string2 = ['a_mean = ' char(string(X1.MXdi)) '\n'];

    % exclude self-connections
    X0.X = C.pdata.X0_f; X0.X(:,di_ind) = 0;
    X1.X = C.pdata.X1_f; X1.X(:,di_ind) = 0;

    % all-inhibitory connections
    X0.I_ind = all(X0.X<0);              X1.I_ind = all(X1.X<0);
    X0.N_I   = length(find(X0.I_ind));   X1.N_I   = length(find(X1.I_ind));
    X0.p_I   = X0.N_I/C.pdata.N_v;       X1.p_I   = X1.N_I/C.pdata.N_v;
    X0.XI    = C.pdata.X0_f(:,X0.I_ind); X1.XI    = C.pdata.X1_f(:,X1.I_ind);
    X0.MXI   = mean(X0.XI,'all');        X1.MXI   = mean(X1.XI,'all');

    I_string1 = 'N_I = 0\n';
    I_string2 = '';
    I_string3 = '';
    if X0.N_I
        I_string1 = ['N_I = '       char(string(X0.N_I)) '\n'];
        I_string2 = ['percent_I = ' char(string(X0.p_I)) '\n'];
        I_string3 = ['a_mean = '    char(string(X0.MXI)) '\n'];
    end

    I_string4 = 'N_I = 0\n';
    I_string5 = '';
    I_string6 = '';
    if X1.N_I
        I_string4 = ['N_I = '       char(string(X1.N_I)) '\n'];
        I_string5 = ['percent_I = ' char(string(X1.p_I)) '\n'];
        I_string6 = ['a_mean = '    char(string(X1.MXI)) '\n'];
    end

    % all-excitatory connections
    X0.E_ind = all(X0.X>0);              X1.E_ind = all(X1.X>0);
    X0.N_E   = length(find(X0.E_ind));   X1.N_E   = length(find(X1.E_ind));
    X0.p_E   = X0.N_E/C.pdata.N_v;       X1.p_E   = X1.N_E/C.pdata.N_v;
    X0.XE    = C.pdata.X0_f(:,X0.E_ind); X1.XE    = C.pdata.X1_f(:,X1.E_ind);
    X0.MXE   = mean(X0.XE,'all');        X1.MXE   = mean(X1.XE,'all');

    E_string1 = 'N_E = 0\n';
    E_string2 = '';
    E_string3 = '';
    if X0.N_E
        E_string1 = ['N_E = '       char(string(X0.N_E)) '\n'];
        E_string2 = ['percent_E = ' char(string(X0.p_E)) '\n'];
        E_string3 = ['a_mean = '    char(string(X0.MXE)) '\n'];
    end

    E_string4 = 'N_E = 0\n';
    E_string5 = '';
    E_string6 = '';
    if X1.N_E
        E_string4 = ['N_E = '       char(string(X1.N_E)) '\n'];
        E_string5 = ['percent_E = ' char(string(X1.p_E)) '\n'];
        E_string6 = ['a_mean = '    char(string(X1.MXE)) '\n'];
    end

    %% Print Info

    outputstr = [ 'E-I Info about ' inputname(1) '\n' ...
                  'SELF-CONNECTIONS\n' ...
                  'Condition 0:\n' di_string1 ...
                  'Condition 1:\n' di_string2 ...
                  'ALL-I-CONNECTIONS\n' ...
                  'Condition 0:\n' ...
                  I_string1 I_string2 I_string3 ...
                  'Condition 1:\n' ...
                  I_string4 I_string5 I_string6 ...
                  'ALL-E-CONNECTIONS\n' ...
                  'Condition 0:\n' ...
                  E_string1 E_string2 E_string3 ...
                  'Condition 1:\n' ...
                  E_string4 E_string5 E_string6 ...
                ];

    display(sprintf(outputstr))
end

end