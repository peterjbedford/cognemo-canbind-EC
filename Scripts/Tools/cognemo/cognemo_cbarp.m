function cognemo_cbarp(C,xset)
%% Preamble
%{
Plots first set of bar plots ('most different' connections by tstat)
%}
%% PREP full

data.IMP_m   = mean(C.class.ccout.IMP,1);
data.IMP_std = std(C.class.ccout.IMP,1);

[data.IMP_m_top,data.IMP_ind] = maxk(data.IMP_m,xset.barp.N_top);
data.IMP_std_top = data.IMP_std(data.IMP_ind);

%% PREP ex

data.ex.IMP_m = data.IMP_m; data.ex.IMP_m(xset.area_ex.c_ind) = 0;

[data.ex.IMP_m_top,data.ex.IMP_ind] = maxk(data.ex.IMP_m,xset.barp.N_top);
data.ex.IMP_std_top = data.IMP_std(data.ex.IMP_ind);

%% Options

options       = xset.barp.options;
options.color = C.xset.pcolor;
options.shared_ind = [];

%% PLOT

cdir = 'Classification bar plots';
if not(isfolder('Results')); mkdir('Results'); end; cd 'Results'
if not(isfolder(cdir)); mkdir(cdir); end; cd(cdir)
if not(isfolder(C.inputname)); mkdir(C.inputname); end; cd(C.inputname)

figs.IMP = figure('Visible','off');
figs.IMP.Units = "points";
figs.IMP.Position = [50 50 xset.barp.sx xset.barp.sy];
cognemo_visualize_feats(data.IMP_m_top,...
                        data.IMP_std_top,...
                        data.IMP_ind,...
                        options)
ax_C = gca;
set(ax_C,'FontSize',9,'TickDir','out')
xlabel(xset.barp.featlabel)
saveas(figs.IMP,[C.inputname '_imp.png'])

%% Plot ex
data.ex.savename = [C.inputname '_imp_ex'];

for i = 1:length(xset.area_ex.label)
    data.ex.savename = [data.ex.savename char(string(xset.area_ex.label(i)))];
    if i < length(xset.area_ex.label)
        data.ex.savename = [data.ex.savename '+'];
    end
end
data.ex.savename = [data.ex.savename '.png'];

figs.ex.IMP = figure('Visible','off');
figs.ex.IMP.Units = "points";
figs.ex.IMP.Position = [50 50 xset.barp.sx xset.barp.sy];
cognemo_visualize_feats(data.ex.IMP_m_top,...
                        data.ex.IMP_std_top,...
                        data.ex.IMP_ind,...
                        options)
ax_C = gca;
set(ax_C,'FontSize',9,'TickDir','out')
xlabel(xset.barp.featlabel)
saveas(figs.ex.IMP,data.ex.savename)

cd ..\..\..

end
