function C = cognemo_mass(C,xset)
%% Preamble
%{
This function wraps tstatistic-based analyses of the connectivity data in
C so that the overall pipeline looks cleaner.
---------------------------------------------------------------------------
INPUTS
---------------------------------------------------------------------------
C:=     structure containing prepared connectivity data pdata.X, which
        split by condition into pdata.X0 and pdata.X1 is fed into a
        t-statistic analysis
xset:=  structure containing analysis settings, including xset.toptions for
        the t-statistic analysis
---------------------------------------------------------------------------
OUTPUTS
---------------------------------------------------------------------------
C:=     same structure as above but with added fields:
            C.tstat:=   contains output from t-statistic analysis across
                        conditions
            C.asy:=     contains output from t-statistic analysis across
                        the diagonal of the connectivity matrix (this is 
                        for analysing 'asymmetry' in EC
%}
%% t-statistic analysis across conditions

C.tstat = cognemo_tstat(C.pdata.X0,C.pdata.X1,...
                        C.pdata.N_c,C.pdata.f_ind,...
                        xset.toptions);

%% asymmetry analysis (across upper and lower triangles)

if ~isfield(C,'inputname')
    C.inputname = inputname(1);
end
if contains(C.inputname,"EC")
    C.asy.pdata = cognemo_getasy(C);      % get upper/lower triangle
    % ttest across upper/lower, condition 0
    C.asy.tstat0 = cognemo_tstat(C.asy.pdata.X0u,C.asy.pdata.X0l,...
                                  C.pdata.N_c,C.asy.pdata.f_ind0,...
                                  xset.toptions);
    % ttest across upper/lower, condition 1
    C.asy.tstat1 = cognemo_tstat(C.asy.pdata.X1u,C.asy.pdata.X1l,...
                                  C.pdata.N_c,C.asy.pdata.f_ind1,...
                                  xset.toptions);
    % ttest of interaction of asy and condition
    C.asy.tstatD = cognemo_tstat(C.asy.pdata.X0u-C.asy.pdata.X0l,...
                                  C.asy.pdata.X1u-C.asy.pdata.X1l,...
                                  C.pdata.N_c,C.asy.pdata.f_ind0,...
                                  xset.toptions);
end

end