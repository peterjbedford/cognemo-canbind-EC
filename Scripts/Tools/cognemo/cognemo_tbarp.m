function cognemo_tbarp(C,xset)
%% Preamble
%{
Plots first set of bar plots ('most different' connections by tstat)
%}
%% PREP full

[~,data.ind_M] = ...
maxk(abs(C.tstat.MXD_f),xset.barp.N_top);
data.MXD = C.tstat.MXD_f(data.ind_M);
data.SXD_M = C.tstat.SXD_f(data.ind_M);

[~,data.ind_T] = ...
maxk(abs(C.tstat.TXD_f),xset.barp.N_top);
data.TXD = C.tstat.TXD_f(data.ind_T);
data.SXD_T = C.tstat.SXD_f(data.ind_T);

%% PREP ex

MXD_ex = C.tstat.MXD_f; MXD_ex(xset.area_ex.c_ind) = 0;

[~,data.ex.ind_M] = ...
maxk(abs(MXD_ex),xset.barp.N_top);
data.ex.MXD = MXD_ex(data.ex.ind_M);
data.ex.SXD_M = C.tstat.SXD_f(data.ex.ind_M);

TXD_ex = C.tstat.TXD_f; TXD_ex(xset.area_ex.c_ind) = 0;

[~,data.ex.ind_T] = ...
maxk(abs(TXD_ex),xset.barp.N_top);
data.ex.TXD = TXD_ex(data.ex.ind_T);
data.ex.SXD_T = C.tstat.SXD_f(data.ex.ind_T);

%% OPTIONS

options = xset.barp.options;
options.color = C.xset.pcolor;
options.shared_ind = [];
options.dir = 1;
if C.inputname == "FC"
    options.dir = 0;
end

%% PLOTS

tsdir = 'T-Statistic bar plots';
if not(isfolder('Results')); mkdir('Results'); end; cd 'Results'
if not(isfolder(tsdir)); mkdir(tsdir); end; cd(tsdir)
if not(isfolder(C.inputname)); mkdir(C.inputname); end; cd(C.inputname)

figs.MXD = figure('Visible','off');
figs.MXD.Units = "points";
figs.MXD.Position = [50 50 xset.barp.sx xset.barp.sy];
cognemo_visualize_feats(data.MXD,...
                        data.SXD_M,...
                        data.ind_M,...
                        options)
ax_C = gca;
set(ax_C,'FontSize',9,'TickDir','out')
xlabel(xset.barp.meanlabel)
saveas(figs.MXD,[C.inputname '_mean.png'])

figs.TXD = figure('Visible','off');
figs.TXD.Units = "points";
figs.TXD.Position = [50 50 xset.barp.sx xset.barp.sy];
cognemo_visualize_feats(data.TXD,...
                        data.SXD_T,...
                        data.ind_T,...
                        options)
ax_C = gca;
set(ax_C,'FontSize',9,'TickDir','out')
xlabel(xset.barp.tstatlabel)
saveas(figs.TXD,[C.inputname '_tstat.png'])

%% PLOT ex

data.ex.MXDsavename = [C.inputname '_mean_ex'];
data.ex.TXDsavename = [C.inputname '_tstat_ex'];

for i = 1:length(xset.area_ex.label)
    data.ex.MXDsavename = [data.ex.MXDsavename char(string(xset.area_ex.label(i)))];
    data.ex.TXDsavename = [data.ex.TXDsavename char(string(xset.area_ex.label(i)))];
    if i < length(xset.area_ex.label)
        data.ex.MXDsavename = [data.ex.MXDsavename '+'];
        data.ex.TXDsavename = [data.ex.TXDsavename '+'];
    end
end
data.ex.MXDsavename = [data.ex.MXDsavename '.png'];
data.ex.TXDsavename = [data.ex.TXDsavename '.png'];

figs.ex.MXD = figure('Visible','off');
figs.ex.MXD.Units = "points";
figs.ex.MXD.Position = [50 50 xset.barp.sx xset.barp.sy];
cognemo_visualize_feats(data.ex.MXD,...
                        data.ex.SXD_M,...
                        data.ex.ind_M,...
                        options)
ax_C = gca;
set(ax_C,'FontSize',9,'TickDir','out')
xlabel(xset.barp.meanlabel)
saveas(figs.ex.MXD,data.ex.MXDsavename)

figs.ex.TXD = figure('Visible','off');
figs.ex.TXD.Units = "points";
figs.ex.TXD.Position = [50 50 xset.barp.sx xset.barp.sy];
cognemo_visualize_feats(data.ex.TXD,...
                        data.ex.SXD_T,...
                        data.ex.ind_T,...
                        options)
ax_C = gca;
set(ax_C,'FontSize',9,'TickDir','out')
xlabel(xset.barp.tstatlabel)
saveas(figs.ex.TXD,data.ex.TXDsavename)

cd ..\..\..

end
