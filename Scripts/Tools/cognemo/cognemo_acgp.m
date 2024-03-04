function cognemo_acgp(C,xset)
%% Preamble
%{
Does the asymmetry connectogram plots
%}
%% SORT

% Condition 0
[CMX0u,~] = cognemo_topk(C.asy.pdata.MX0u_f,xset.cgp.N_top,1);
data.MX0u = ...
    cognemo_cg_prep(cognemo_symmtx(CMX0u,1),xset.cgp.ind); clear CMX0u
[CMX0l,~] = cognemo_topk(C.asy.pdata.MX0l_f,xset.cgp.N_top,1);
data.MX0l = ...
    cognemo_cg_prep(cognemo_symmtx(CMX0l,1),xset.cgp.ind); clear CMX0l
[CTX0D,~] = cognemo_topk(C.asy.tstat0.TXD_f,xset.cgp.N_top,1);
data.TX0D = ...
    cognemo_cg_prep(cognemo_symmtx(CTX0D,1),xset.cgp.ind); clear CTX0D

% Condition 1
[CMX1u,~] = cognemo_topk(C.asy.pdata.MX1u_f,xset.cgp.N_top,1);
data.MX1u = ...
    cognemo_cg_prep(cognemo_symmtx(CMX1u,1),xset.cgp.ind); clear CMX1u
[CMX1l,~] = cognemo_topk(C.asy.pdata.MX1l_f,xset.cgp.N_top,1);
data.MX1l = ...
    cognemo_cg_prep(cognemo_symmtx(CMX1l,1),xset.cgp.ind); clear CMX1l
[CTX1D,~] = cognemo_topk(C.asy.tstat1.TXD_f,xset.cgp.N_top,1);
data.TX1D = ...
    cognemo_cg_prep(cognemo_symmtx(CTX1D,1),xset.cgp.ind); clear CTX1D


%% PLOT

asydir = 'Asymmetry connectogram plots';
if not(isfolder('Results')); mkdir('Results'); end; cd 'Results'
if not(isfolder(asydir)); mkdir(asydir); end; cd(asydir)
if not(isfolder(C.inputname)); mkdir(C.inputname); end; cd(C.inputname)

% Condition 0
figs.MX0u = figure('Visible','off');
figs.MX0u.Position = [50 50 xset.cgp.sx xset.cgp.sy];
circularGraph(data.MX0u,'Label',cellstr(xset.cgp.label));
saveas(figs.MX0u,[C.inputname '_cond0_triu.png']) 
figs.MX0l = figure('Visible','off');
figs.MX0l.Position = [50 50 xset.cgp.sx xset.cgp.sy];
circularGraph(data.MX0l,'Label',cellstr(xset.cgp.label));
saveas(figs.MX0l,[C.inputname '_cond0_tril.png'])
figs.TX0D = figure('Visible','off');
figs.TX0D.Position = [50 50 xset.cgp.sx xset.cgp.sy];
circularGraph(data.TX0D,'Label',cellstr(xset.cgp.label));
saveas(figs.TX0D,[C.inputname '_cond0_triu-tril.png'])

% Condition 0
figs.MX1u = figure('Visible','off');
figs.MX1u.Position = [50 50 xset.cgp.sx xset.cgp.sy];
circularGraph(data.MX1u,'Label',cellstr(xset.cgp.label));
saveas(figs.MX1u,[C.inputname '_cond1_triu.png'])
figs.MX1l = figure('Visible','off');
figs.MX1l.Position = [50 50 xset.cgp.sx xset.cgp.sy];
circularGraph(data.MX1l,'Label',cellstr(xset.cgp.label));
saveas(figs.MX1l,[C.inputname '_cond1_tril.png'])
figs.TX1D = figure('Visible','off');
figs.TX1D.Position = [50 50 xset.cgp.sx xset.cgp.sy];
circularGraph(data.TX1D,'Label',cellstr(xset.cgp.label));
saveas(figs.TX1D,[C.inputname '_cond1_triu-tril.png'])

cd ..\..\..

end