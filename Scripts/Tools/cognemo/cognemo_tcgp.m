function cognemo_tcgp(C,xset)
%% Preamble
%{
Does the first set of connectogram plots:
- mean of condition 1
- mean of condition 2
- tstat of difference: condition 1 - condition 2
%}
%% SORT
method = C.xset.symmethod;

[CMX0,~] = cognemo_topk(C.pdata.MX0_f,xset.cgp.N_top,1);
data.MX0 = ...
    cognemo_cg_prep(cognemo_symmtx(CMX0,method),xset.cgp.ind); clear CMX0
[CMX1,~] = cognemo_topk(C.pdata.MX1_f,xset.cgp.N_top,1);
data.MX1 = ...
    cognemo_cg_prep(cognemo_symmtx(CMX1,method),xset.cgp.ind); clear CMX1
[CTXD,~] = cognemo_topk(C.tstat.TXD_f,xset.cgp.N_top,1);
data.TXD = ...
    cognemo_cg_prep(cognemo_symmtx(CTXD,method),xset.cgp.ind); clear CTXD
%% PLOT

tsdir = 'T-Statistic connectogram plots';
if not(isfolder('Results')); mkdir('Results'); end; cd 'Results'
if not(isfolder(tsdir)); mkdir(tsdir); end; cd(tsdir)
if not(isfolder(C.inputname)); mkdir(C.inputname); end; cd(C.inputname)

figs.MX0 = figure('Visible','off');
figs.MX0.Position = [50 50 xset.cgp.sx xset.cgp.sy];
circularGraph(data.MX0,'Label',cellstr(xset.cgp.label));
figs.MX0_T = cognemo_ringlabel;
saveas(figs.MX0,[C.inputname '_LSD_mean.png'])
figs.MX1 = figure('Visible','off');
figs.MX1.Position = [50 50 xset.cgp.sx xset.cgp.sy];
circularGraph(data.MX1,'Label',cellstr(xset.cgp.label));
figs.MX1_T = cognemo_ringlabel;
saveas(figs.MX1,[C.inputname '_PLA_mean.png'])
figs.TXD = figure('Visible','off');
figs.TXD.Position = [50 50 xset.cgp.sx xset.cgp.sy];
circularGraph(data.TXD,'Label',cellstr(xset.cgp.label));
figs.TXD_T = cognemo_ringlabel;
saveas(figs.TXD,[C.inputname '_LSD-PLA_mean.png'])

cd ..\..\..

end