%% CANBIND ANALYSIS--MDD vs HC

% Step 2: Performance analysis

% NOTE: Control = 0 | MDD = 1

%% Add tools to path

addpath('Tools/circularGraph')
addpath('Tools/cognemo')
addpath('Tools/cognemoPL')
addpath('Tools/subaxis')

%% Save images? Display images?

pset.disp = 'on';
pset.save = 1;
pset.atlas = "AAL";

% connection label options
pset.coptions.dir = 1;
pset.coptions.exstr = "_";

%% Import RegionLabels table

labeltab = readtable('RegionLabels_Nov22.xlsx','Sheet','Sheet1','Range','A1:M271');
%}
%% Convert RegionLabels variables to matlab-friendly and save
%{
rlabel.ind_sort = table2array(labeltab(:,1));
rlabel.AALlabel_sort = string(table2array(labeltab(:,3)));
rlabel.BRlabel_sort = string(table2array(labeltab(:,4)));
rlabel.fnlabel_sort = string(table2array(labeltab(:,5)));
rlabel.fnlabsh_sort = string(table2array(labeltab(:,6)));
rlabel.color_sort  = string(table2array(labeltab(:,7)));
rlabel.rlabel_sort = rlabel.AALlabel_sort;
%}
rlabel.fn_unsort = table2array(labeltab(:,2));
[~,rlabel.fn_sort] = sort(rlabel.fn_unsort,'ascend');
rlabel.ar_unsort = table2array(labeltab(:,3));
[~,rlabel.ar_sort] = sort(rlabel.ar_unsort,'ascend');
rlabel.AAL.og  = string(table2array(labeltab(:,5)));
rlabel.bdm.og  = string(table2array(labeltab(:,6)));
rlabel.fnl.og  = string(table2array(labeltab(:,7)));
rlabel.fns.og  = string(table2array(labeltab(:,8)));
rlabel.arl.og  = string(table2array(labeltab(:,9)));
rlabel.ars.og  = string(table2array(labeltab(:,10)));
rlabel.fnc.og  = string(table2array(labeltab(:,11)));
rlabel.arc.og  = string(table2array(labeltab(:,12)));

N_r = size(rlabel.fnc.og,1);
fnc2 = zeros([N_r,3]); arc2 = fnc2;
for i = 1:N_r
    fnc_i = erase(rlabel.fnc.og(i),["[","]"]);
    arc_i = erase(rlabel.arc.og(i),["[","]"]);
    fnc_i = str2double(strsplit(string(fnc_i)));
    arc_i = str2double(strsplit(string(arc_i)));
    fnc2(i,:) = fnc_i;
    arc2(i,:) = arc_i;
end
rlabel = rmfield(rlabel,{'fnc','arc'});
rlabel.fnc.og = fnc2; rlabel.arc.og = arc2; clear i fnc2 fnc_i arc2 arc_i N_r

save('rlabel.mat','rlabel')
clear labeltab rlabel

%% Import rlabel variables
load('rlabel.mat');
if pset.atlas == "AAL"; rlabel.rlabel.og = rlabel.AAL.og;
else rlabel.rlabel.og = rlabel.bdm.og; end
save('rlabel.mat','rlabel')
%% Create Connection Labels and save
%{
N_r = length(rlabel.rlabel_og); N_c = N_r*N_r;
options.dir = 1; options.exstr = "_"; options.areacolors = rlabel.fnc_og;
clabel.clabel = cognemo_ind2clabel(1:N_c,rlabel.rlabel_og',options);
clabel.clabel_sort = reshape(clabel.clabel,[N_r,N_r]);
clabel.clabel_sort = clabel.clabel_sort(rlabel.ind_sort,rlabel.ind_sort);
clabel.clabel_sort = reshape(clabel.clabel_sort,[1,N_c]);
clear N_r N_c 
save('clabel.mat','clabel'); clear clabel ind_unsort
%}
N_r = length(rlabel.rlabel.og); N_c = N_r*N_r;
clabel.clabel.og = cognemo_ind2clabel(1:N_c,rlabel.rlabel.og',pset.coptions);
clabelsq = reshape(clabel.clabel.og,[N_r,N_r]);
clabel.clabel.fn = clabelsq(rlabel.fn_sort,rlabel.fn_sort);
clabel.clabel.fn = reshape(clabel.clabel.fn,[1,N_c]);
clear N_r N_c
save('clabel.mat','clabel'); clear clabel clabelsq

%% Import clabel variables
load('clabel.mat')
%% Figure settings

pset.pcolor = [255 181 0]/255;
pset.units = 'inches';
pset.dispconv = 0.48;
pset.tstatlabel = 't-statistic of MDD-HC';
pset.tstatsclabel = 't-statistic of MDD-HC (SC)';
pset.featlabel = 'Most Discriminative';
pset.featsclabel = 'Most Discriminative (SC)';

% hyperparameter plot settings
pset.hyset.figpos = [0 0 8 8]*pset.dispconv;
pset.hyset.fontsize = 12;
pset.hyset.linewidth = 2;
N_color = 1000; pset.hyset.colormap = zeros([N_color,3]); 
pset.hyset.colormap(:,1) = (1:N_color)*pset.pcolor(1)/N_color;
pset.hyset.colormap(:,2) = (1:N_color)*pset.pcolor(2)/N_color;
pset.hyset.colormap(:,3) = (1:N_color)*pset.pcolor(3)/N_color; clear N_color

% connectogram settings
pset.cgset.N_top = 250;
pset.cgset.figpos = [0 0 18 18]*pset.dispconv;
pset.cgset.rlabel = rlabel; N_r = length(pset.cgset.rlabel.ar_sort);
indf25 = 1:ceil(0.25*N_r); indl75 = (ceil(0.25*N_r)+1):N_r;
indf75 = 1:length(indl75); indl25 = (length(indf75)+1):N_r;
pset.cgset.rlabel.ar_sort(indl25) = rlabel.ar_sort(indf25);
pset.cgset.rlabel.ar_sort(indf75) = rlabel.ar_sort(indl75);
clear indf25 indl75 indf75 indl25
pset.cgset.fn.colormap = rlabel.fnc.og(rlabel.fn_sort,:)./255;
pset.cgset.fn.ringoptions.radius = 1.025;
pset.cgset.fn.ringoptions.fontsize = 10;
pset.cgset.fn.ringoptions.rotate = 1;
pset.cgset.ar.colormap = rlabel.arc.og(pset.cgset.rlabel.ar_sort,:)./255;
pset.cgset.ar.ringoptions.radius = 1.025;
pset.cgset.ar.ringoptions.fontsize = 10;
pset.cgset.ar.ringoptions.rotate = 0;
pset.cgset.labelblank = cellstr(repmat(" ",[270,1]));
pset.cgset.dim = [0.05 0.125 0.4 0.775];

% bar plot (many) settings
pset.bmset.N_top = 50;
pset.bmset.figpos = [0 0 11 9]*pset.dispconv;
pset.bmset.fontsize = 8;
options.arealabel = rlabel.ars.og;
options.areacolors = rlabel.fnc.og;
options.pcolor = pset.pcolor;
options.dir = 1;
pset.bmset.options = options; clear options
legoptions.arealabels = rlabel.fns.og(rlabel.fn_sort);
legoptions.areacolor = rlabel.fnc.og(rlabel.fn_sort);
legoptions.fontsize = 8;
legoptions.linespacing = 1.5;
pset.bmset.legoptions = legoptions; clear legoptions

% bar plot (few) settings
pset.bfset.N_top = 10;
pset.bfset.figpos = [0 0 8 5]*pset.dispconv;
pset.bfset.fontsize = 8;
options.N_top = pset.bfset.N_top;
options.clabel = clabel.clabel;
options.clabeldir = clabel.clabel;
options.shared_ind = [];
options.shared_color = [];
options.errbar = 'std';
options.dir = 1;
pset.bfset.options = options; clear options

% network method settings
pset.nmset = pset.hyset;
pset.nmset.figpos = [0 0 12 8]*pset.dispconv;

clear rlabel clabel

%% Input data

load('in.mat')
EC.in.X = X; EC.in.T = T; EC.in.X0 = X(~logical(T),:); EC.in.X1 = X(logical(T),:);
[EC.in.N_o,EC.in.N_c] = size(EC.in.X); EC.in.N_r = sqrt(EC.in.N_c);
EC.in.ind_sc = find(logical(reshape(eye(EC.in.N_r),[1,EC.in.N_c])));
EC.in.ind_nsc = find(~logical(reshape(eye(EC.in.N_r),[1,EC.in.N_c])));

%% BASIC STATS PLOTS

% basic stats
EC.stats.MX0 = mean(EC.in.X0,1); EC.stats.SX0 = std(EC.in.X0,1);
EC.stats.MX0_sc  = EC.stats.MX0; EC.stats.MX0_sc(EC.in.ind_nsc) = 0;
EC.stats.MX0_nsc = EC.stats.MX0; EC.stats.MX0_nsc(EC.in.ind_sc) = 0;
[~,EC.stats.top_ind_MX0] = sort(abs(EC.stats.MX0),'descend');
[~,EC.stats.top_ind_MX0_sc] = sort(abs(EC.stats.MX0_sc),'descend');
[~,EC.stats.top_ind_MX0_nsc] = sort(abs(EC.stats.MX0_nsc),'descend');
EC.stats.MX1 = mean(EC.in.X1,1); EC.stats.SX1 = std(EC.in.X1,1);
EC.stats.MX1_sc  = EC.stats.MX1; EC.stats.MX1_sc(EC.in.ind_nsc) = 0;
EC.stats.MX1_nsc = EC.stats.MX1; EC.stats.MX1_nsc(EC.in.ind_sc) = 0;
[~,EC.stats.top_ind_MX1] = sort(abs(EC.stats.MX1),'descend');
[~,EC.stats.top_ind_MX1_sc] = sort(abs(EC.stats.MX1_sc),'descend');
[~,EC.stats.top_ind_MX1_nsc] = sort(abs(EC.stats.MX1_nsc),'descend');
% ttest; MDD - control
[h,~,~,stats] = ttest2(EC.in.X1,EC.in.X0); EC.stats.N_sig = length(find(h));
EC.stats.TXD = stats.tstat; EC.stats.SXD = stats.sd; clear stats h
EC.stats.TXD_sc  = EC.stats.TXD; EC.stats.TXD_sc(EC.in.ind_nsc) = 0;
EC.stats.TXD_nsc = EC.stats.TXD; EC.stats.TXD_nsc(EC.in.ind_sc) = 0;
[~,EC.stats.top_ind_TXD] = sort(abs(EC.stats.TXD),'descend');
[~,EC.stats.top_ind_TXD_sc] = sort(abs(EC.stats.TXD_sc),'descend');
[~,EC.stats.top_ind_TXD_nsc] = sort(abs(EC.stats.TXD_nsc),'descend');

%% Figure 1: Settings

fig1.set.pos = [0 0 1500 1500];
fig1.set.N_top = 250;
fig1.set.cgp.labels = cellstr(repmat([" "],[1,N_r]));
fig1.set.cbar.fontsize = 16;
fig1.set.cbar.colormin = [0 201 239]./255;
fig1.set.cbar.colormax = [255 181 0]./255;
fig1.set.title.labels{1,1} = '\bfA. \rm mean MDD';
fig1.set.title.labels{1,2} = '\bfB. \rm mean HC';
fig1.set.title.labels{1,3} = '\bfC. \rm t-value MDD-HC';
fig1.set.title.labels{2,1} = '\bfD. \rm mean MDD';
fig1.set.title.labels{2,2} = '\bfE. \rm mean HC';
fig1.set.title.labels{2,3} = '\bfF. \rm t-value MDD-HC';
fig1.set.title.fontsize = 16;

%% Figure 1: Data
N_c = N_r*N_r;
di_ind = logical(reshape(eye(N_r),[1,N_c]));

fig1.data.in = cell(2,3);
fig1.data.in{1,1} = EC.stats.MX0; fig1.data.in{2,1} = EC.stats.MX0;
fig1.data.in{1,2} = EC.stats.MX1; fig1.data.in{2,2} = EC.stats.MX1;
fig1.data.in{1,3} = EC.stats.TXD; fig1.data.in{2,3} = EC.stats.TXD;
fig.data.out = cell(2,3);
ymins = zeros([1,6]); ymaxs = zeros([1,6]);

for j = 1:6
    if j<=3
        i_data = 1;
        j_data = j;
        r_sort = pset.cgset.rlabel.fn_sort;
    else
        i_data = 2;
        j_data = j - 3;
        r_sort = pset.cgset.rlabel.ar_sort;
    end
    
    datatemp = fig1.data.in{i_data,j_data};
    % if...
    datatemp(:,di_ind) = 0;
    [datatemp,~] = cognemo_topk(datatemp,fig1.set.N_top,1);
    fig1.data.out{i_data,j_data} = ...
        cognemo_cg_prep(cognemo_symmtx(datatemp,2),r_sort); clear datatemp
    datatemp = fig1.data.out{i_data,j_data};
    ymins(j) = min(datatemp(logical(datatemp)));
    ymaxs(j) = max(datatemp(logical(datatemp)));
    
    clear datatemp ymin ymax i_data j_data r_sort

end; clear j

% subplot groups
fig1.set.i_groupa = [1,2,4,5]; fig1.set.i_groupb = [3,6];
fig1.data.min_groupa = min(ymins(fig1.set.i_groupa));
fig1.data.min_groupb = min(ymins(fig1.set.i_groupb));
fig1.data.max_groupa = max(ymaxs(fig1.set.i_groupa));
fig1.data.max_groupb = max(ymaxs(fig1.set.i_groupb));
clear ymins ymaxs

%% Figure 1: Plot
cmd.show = 1;
fig1.set.fprop.Name = {'OuterPosition','Color','Visible'};
fig1.set.fprop.Value = {fig1.set.pos,'k',cmd.show};
fig1.set.aprop.Name = {'Color','XColor','YColor','LineWidth'};
fig1.set.aprop.Value = {'k','w','w',2};
fig1.draw.fig = figure();
set(fig1.draw.fig,fig1.set.fprop.Name,fig1.set.fprop.Value);
adjust = 0.01;

for j = 1:6
    % top row
    if j<=3
        i_data = 1;
        j_data = j;
        f1_colors = pset.cgset.fn.colormap;
        f1_labels = pset.cgset.rlabel.fns.og(pset.cgset.rlabel.fn_sort);
        f1_ringoptions = pset.cgset.fn.ringoptions;
        axpos = [ 0.33*(j_data-1) 0.5 0.33 0.5 ];
        cbarpos = [ 0.33*(j_data-1)+adjust 0.5+adjust 0.04 0.125 ];
        titlepos = [ 0.33*(j_data-1)+adjust 0.5 0.33 0.5-adjust ];
        ringfun = @cognemo_ringlabel;
    % bottom row
    else
        i_data = 2;
        j_data = j - 3;
        f1_colors = pset.cgset.ar.colormap;
        f1_labels = pset.cgset.rlabel.ars.og(pset.cgset.rlabel.ar_sort);
        f1_ringoptions = pset.cgset.ar.ringoptions;
        axpos = [ 0.33*(j_data-1) 0 0.33 0.5 ];
        cbarpos = [ 0.33*(j_data-1)+adjust 0+adjust 0.04 0.125 ];
        titlepos = [ 0.33*(j_data-1)+adjust 0 0.33 0.5-adjust ];
        ringfun = @cognemo_ringlabel2;
    end
    % A,B,D,E have same scale
    if ismember(j,fig1.set.i_groupa)
        f1_lims = [fig1.data.min_groupa,fig1.data.max_groupa];
    % C,F have same scale
    else
        f1_lims = [fig1.data.min_groupb,fig1.data.max_groupb];
    end
    
    axes('OuterPosition',axpos,'Color','k');
    % Connectogram
    circularGraph_lsd(fig1.data.out{i_data,j_data},...
        'Label',fig1.set.cgp.labels,...
        'Areacolors',f1_colors,...
        'Limits',f1_lims);
    ringfun(f1_labels,f1_ringoptions);
    % Colorbar
    cognemo_alphabar(f1_lims,cbarpos,fig1.set.cbar);
    annotation('textbox',titlepos,...
               'String',fig1.set.title.labels{i_data,j_data},...
               'Margin',0,...
               'Color','w','FontSize',fig1.set.title.fontsize,...
               'LineStyle','none');
    
    
    clear i_data j_data f1_colors f1_labels f1_ringoptions f1_lims

end

clear adjust

%% cgp -- network views
%% cgp meanHC (MX0) prep
cg.meanHC.data = zeros(size(EC.stats.MX0));
cg.meanHC.data(EC.stats.top_ind_MX0_nsc(1:pset.cgset.N_top)) = ...
    EC.stats.MX0_nsc(EC.stats.top_ind_MX0_nsc(1:pset.cgset.N_top));
cg.meanHC.data = cognemo_symmtx(cg.meanHC.data,2);
cg.meanHC.data = cognemo_cg_prep(cg.meanHC.data,pset.cgset.rlabel.fn_sort);
%% cgp meanHC (MX0) output plot
cg.meanHC.fig = figure('Units',pset.units,...
                       'Visible',pset.disp,'Position',pset.cgset.figpos);
circularGraph(cg.meanHC.data,'Label',pset.cgset.labelblank,'Areacolors',pset.cgset.fn.colormap);
cognemo_ringlabel(pset.cgset.rlabel.fns.og(pset.cgset.rlabel.fn_sort),pset.cgset.fn.ringoptions)
if pset.save; cg.meanHC.fname = 'cgfn_meanHC.png'; cd ..; cd Figures
saveas(cg.meanHC.fig,cg.meanHC.fname); cd ..; cd Pipeline; end
%% cgp meanMDD (MX1) prep
cg.meanMDD.data = zeros(size(EC.stats.MX1));
cg.meanMDD.data(EC.stats.top_ind_MX1_nsc(1:pset.cgset.N_top)) = ...
    EC.stats.MX1_nsc(EC.stats.top_ind_MX1_nsc(1:pset.cgset.N_top));
cg.meanMDD.data = cognemo_symmtx(cg.meanMDD.data,2);
cg.meanMDD.data = cognemo_cg_prep(cg.meanMDD.data,pset.cgset.rlabel.fn_sort);
%% cgp meanHC (MX1) output plot
cg.meanMDD.fig = figure('Units',pset.units,...
                       'Visible',pset.disp,'Position',pset.cgset.figpos);
circularGraph(cg.meanMDD.data,'Label',pset.cgset.labelblank,'Areacolors',pset.cgset.fn.colormap);
cognemo_ringlabel(pset.cgset.rlabel.fns.og(pset.cgset.rlabel.fn_sort),pset.cgset.fn.ringoptions)
if pset.save; cg.meanMDD.fname = 'cgfn_meanMDD.png'; cd ..; cd Figures
saveas(cg.meanMDD.fig,cg.meanMDD.fname); cd ..; cd Pipeline; end

%% cgp tstat prep
cg.tstat.data = zeros(size(EC.stats.TXD));
cg.tstat.data(EC.stats.top_ind_TXD_nsc(1:pset.cgset.N_top)) = ...
    EC.stats.TXD_nsc(EC.stats.top_ind_TXD_nsc(1:pset.cgset.N_top));
cg.tstat.data = cognemo_symmtx(cg.tstat.data,2);
cg.tstat.data = cognemo_cg_prep(cg.tstat.data,pset.cgset.rlabel.fn_sort);
%% cgp tstat output plot
cg.tstat.fig = figure('Units',pset.units,...
                      'Visible',pset.disp,'Position',pset.cgset.figpos);
circularGraph(cg.tstat.data,'Label',pset.cgset.labelblank,'Areacolors',pset.cgset.fn.colormap);
cognemo_ringlabel(pset.cgset.rlabel.fns.og(pset.cgset.rlabel.fn_sort),pset.cgset.fn.ringoptions)
if pset.save; cg.tstat.fname = 'cgfn_tstat_MDD-HC.png'; cd ..; cd Figures
saveas(cg.tstat.fig,cg.tstat.fname); cd ..; cd Pipeline; end

%% cgp -- anatomical area views
pset.cgset.ar.ringoptions.rotate = 0;
%% cgp meanHC (MX0) prep
cg.meanHC.data = zeros(size(EC.stats.MX0));
cg.meanHC.data(EC.stats.top_ind_MX0_nsc(1:pset.cgset.N_top)) = ...
    EC.stats.MX0_nsc(EC.stats.top_ind_MX0_nsc(1:pset.cgset.N_top));
cg.meanHC.data = cognemo_symmtx(cg.meanHC.data,2);
cg.meanHC.data = cognemo_cg_prep(cg.meanHC.data,pset.cgset.rlabel.ar_sort);
%% cgp meanHC (MX0) output plot
cg.meanHC.fig = figure('Units',pset.units,...
                       'Visible',pset.disp,'Position',pset.cgset.figpos);
circularGraph(cg.meanHC.data,'Label',pset.cgset.labelblank,'Areacolors',pset.cgset.ar.colormap);
cognemo_ringlabel2(pset.cgset.rlabel.ars.og(pset.cgset.rlabel.ar_sort),pset.cgset.ar.ringoptions)
if pset.save; cg.meanHC.fname = 'cgar_meanHC.png'; cd ..; cd Figures
saveas(cg.meanHC.fig,cg.meanHC.fname); cd ..; cd Pipeline; end
%% cgp meanMDD (MX1) prep
cg.meanMDD.data = zeros(size(EC.stats.MX1));
cg.meanMDD.data(EC.stats.top_ind_MX1_nsc(1:pset.cgset.N_top)) = ...
    EC.stats.MX1_nsc(EC.stats.top_ind_MX1_nsc(1:pset.cgset.N_top));
cg.meanMDD.data = cognemo_symmtx(cg.meanMDD.data,2);
cg.meanMDD.data = cognemo_cg_prep(cg.meanMDD.data,pset.cgset.rlabel.ar_sort);
%% cgp meanMDD (MX1) output plot
cg.meanMDD.fig = figure('Units',pset.units,...
                       'Visible',pset.disp,'Position',pset.cgset.figpos);
circularGraph(cg.meanMDD.data,'Label',pset.cgset.labelblank,'Areacolors',pset.cgset.ar.colormap);
cognemo_ringlabel2(pset.cgset.rlabel.ars.og(pset.cgset.rlabel.ar_sort),pset.cgset.ar.ringoptions)
if pset.save; cg.meanMDD.fname = 'cgar_meanMDD.png'; cd ..; cd Figures
saveas(cg.meanMDD.fig,cg.meanMDD.fname); cd ..; cd Pipeline; end

%% cgp tstat prep
cg.tstat.data = zeros(size(EC.stats.TXD));
cg.tstat.data(EC.stats.top_ind_TXD_nsc(1:pset.cgset.N_top)) = ...
    EC.stats.TXD_nsc(EC.stats.top_ind_TXD_nsc(1:pset.cgset.N_top));
cg.tstat.data = cognemo_symmtx(cg.tstat.data,2);
cg.tstat.data = cognemo_cg_prep(cg.tstat.data,pset.cgset.rlabel.ar_sort);
%% cgp tstat output plot
cg.tstat.fig = figure('Units',pset.units,...
                      'Visible',pset.disp,'Position',pset.cgset.figpos);
circularGraph(cg.tstat.data,'Label',pset.cgset.labelblank,'Areacolors',pset.cgset.ar.colormap);
cognemo_ringlabel2(pset.cgset.rlabel.ars.og(pset.cgset.rlabel.ar_sort),pset.cgset.ar.ringoptions)
if pset.save; cg.tstat.fname = 'cgar_tstat_MDD-HC.png'; cd ..; cd Figures
saveas(cg.tstat.fig,cg.tstat.fname); cd ..; cd Pipeline; end

%%
%% bfp tstat (TXD) prep
bf.tstat.data.ind = EC.stats.top_ind_TXD(1:pset.bfset.N_top);
bf.tstat.data.TXD = ...
    EC.stats.TXD(bf.tstat.data.ind);
bf.tstat.data.SXD = ...
    EC.stats.SXD(bf.tstat.data.ind);
%% bfp tstat (TXD) plot
bf.tstat.fig = figure('Units',pset.units,...
                      'Visible',pset.disp,'Position',pset.bfset.figpos);
cognemo_visualize_feats(bf.tstat.data.TXD,...
                        bf.tstat.data.SXD,...
                        bf.tstat.data.ind,...
                        pset.bfset.options); clear options
set(gca,'FontSize',pset.bfset.fontsize,'TickDir','in')
title(pset.tstatlabel,'Color','w','FontSize',1.25*get(gca,'FontSize'))
xlabel('\rightarrow MDD > HC')
if pset.save; bf.tstat.fname = 'bf_tstat_MDD-HC.png'; cd ..; cd Figures
saveas(bf.tstat.fig,bf.tstat.fname); cd ..; cd Pipeline; end
%% bfp tstatsc (TXD_sc) prep
bf.tstatsc.data.ind = EC.stats.top_ind_TXD_sc(1:pset.bfset.N_top);
bf.tstatsc.data.TXD = ...
    EC.stats.TXD(bf.tstatsc.data.ind);
bf.tstatsc.data.SXD = ...
    EC.stats.SXD(bf.tstatsc.data.ind);
%% bfp tstatsc (TXD_sc) plot
bf.tstatsc.fig = figure('Units',pset.units,...
                        'Visible',pset.disp,'Position',pset.bfset.figpos);
cognemo_visualize_feats(bf.tstatsc.data.TXD,...
                        bf.tstatsc.data.SXD,...
                        bf.tstatsc.data.ind,...
                        pset.bfset.options); clear options
set(gca,'FontSize',pset.bfset.fontsize,'TickDir','in')
title(pset.tstatsclabel,'Color','w','FontSize',1.25*get(gca,'FontSize'))
xlabel('MDD < HC \leftarrow          \rightarrow MDD > HC')
if pset.save; bf.tstatsc.fname = 'bf_tstatsc_MDD-HC.png'; cd ..; cd Figures
saveas(bf.tstatsc.fig,bf.tstatsc.fname); cd ..; cd Pipeline; end

%% bmp tstat (TXD) prep
bm.tstat.data.ind = EC.stats.top_ind_TXD(1:pset.bmset.N_top);
bm.tstat.data.TXD = ...
    EC.stats.TXD(bm.tstat.data.ind);

%% bmp tstat (TXD) plot
bm.tstat.fig = figure('Units',pset.units,...
                      'Visible',pset.disp,'Position',pset.bmset.figpos,...
                      'inverthardcopy','off','color','k');
cognemo_areabar(bm.tstat.data.TXD,bm.tstat.data.ind,pset.bmset.options);
xlabel('MDD < HC \leftarrow          \rightarrow MDD > HC')
set(gca,'Position',[0.075 0.125 0.5375 0.75],...
        'FontSize',pset.bmset.fontsize);
title(pset.tstatlabel,'Color','w','FontSize',1.25*get(gca,'FontSize'))
% annotation('textbox',[0.0125 0.025 0.6375 0.95],...
%            'String',"",...
%            'Color','w','EdgeColor','w','FontSize',24,'LineWidth',2);
cognemo_areabar_legend([0.65 0.15 0.4 0.75],pset.bmset.legoptions)
if pset.save; bm.tstat.fname = 'bm_tstat_MDD-HC.png'; cd ..; cd Figures
saveas(bm.tstat.fig,bm.tstat.fname); cd ..; cd Pipeline; end

%% bmp tstatsc (TXD_sc) prep
bm.tstatsc.data.ind = EC.stats.top_ind_TXD_sc(1:pset.bmset.N_top);
bm.tstatsc.data.TXD = ...
    EC.stats.TXD(bm.tstatsc.data.ind);

%% bmp tstatsc (TXD_sc) plot
bm.tstatsc.fig = figure('Units',pset.units,...
                        'Visible',pset.disp,'Position',pset.bmset.figpos,...
                        'inverthardcopy','off','color','k');
cognemo_areabar(bm.tstatsc.data.TXD,bm.tstatsc.data.ind,pset.bmset.options);
xlabel('MDD < HC \leftarrow          \rightarrow MDD > HC')
set(gca,'Position',[0.075 0.125 0.5375 0.75],...
        'FontSize',pset.bmset.fontsize);
title(pset.tstatsclabel,'Color','w','FontSize',1.25*get(gca,'FontSize'))
% annotation('textbox',[0.0125 0.025 0.6375 0.95],...
%            'String',"",...
%            'Color','w','EdgeColor','w','FontSize',24,'LineWidth',2);
cognemo_areabar_legend([0.65 0.15 0.4 0.75],pset.bmset.legoptions)
if pset.save; bm.tstatsc.fname = 'bm_tstatsc_MDD-HC.png'; cd ..; cd Figures
saveas(bm.tstatsc.fig,bm.tstatsc.fname); cd ..; cd Pipeline; end


%% Open output data

cd("NCVRF_18-Dec-2021 16-10-29")
load('NCVRF_18-Dec-2021 16-10-29.mat')
cd ..

%% HYPERPARAMETER PLOT

%% Data

[~,I_N_f]   = sort(out.O.N_f,'ascend');
[~,I_R_ttf] = sort(out.O.R_ttf,'ascend');

OCgrid = mean(out.I.OC,4); OCgrid = mean(OCgrid,3);

hy.Xdata = out.O.N_f;
hy.Ydata = out.O.R_ttf;
hy.Cdata = OCgrid;

%% Apply sort
hy.Xdata = hy.Xdata(I_N_f);
hy.Ydata = hy.Ydata(I_R_ttf);
hy.Cdata = hy.Cdata(I_R_ttf,I_N_f); Cdata_norm = normalize(hy.Cdata,'range');
Cdata_temp = zeros([size(hy.Cdata,1),size(hy.Cdata,2),3]);
Cdata_temp(:,:,1) = Cdata_norm*pset.pcolor(1);
Cdata_temp(:,:,2) = Cdata_norm*pset.pcolor(2);
Cdata_temp(:,:,3) = Cdata_norm*pset.pcolor(3);
% hy.Cdata = Cdata_temp;
clear Cdata_temp C_data_norm


%% Optimization grid plot - figure

% Labels
hy.xlabel = "Number of features";
hy.ylabel = "Ratio of trees-to-features";

% Figure
hy.fig = figure('Units',pset.units,...
                'Visible',pset.disp,...
                'OuterPosition',pset.hyset.figpos);
imagesc(hy.Cdata); colormap(pset.hyset.colormap);
colorbar('Color','w','LineWidth',2);
hy.xticks = 1:length(get(gca,'XTick'));
hy.yticks = hy.xticks;
xticks(hy.xticks); xticklabels(string(hy.Xdata)); xlabel(hy.xlabel);
yticks(hy.yticks); yticklabels(string(hy.Ydata)); ylabel(hy.ylabel);
% Style
set(gcf,'Color','k','InvertHardCopy','off');
set(gca,'XColor','w','YColor','w',...
        'FontSize',pset.hyset.fontsize,...
        'LineWidth',pset.hyset.linewidth)
if pset.save; hy.fname = 'hyperparameters.png'; cd ..; cd Figures
saveas(hy.fig,hy.fname); cd ..; cd Pipeline; end


%% Permutation test


%% load outer fold models
cd("NCVRF_18-Dec-2021 16-10-29/outermdl")

%%
fnames = strsplit(ls('*.mat'));
fnames = sort(fnames);
fnames = fnames(2:length(fnames));

% permuting T
N_perm = 1000; rng(1234);
T0     = [ ones([1,174/2]) zeros([1,174/2]) ];
T_all   = zeros([N_perm,174]); OC_all = zeros([N_perm,options.KO]);
for n = 1:N_perm
    I_n = randperm(size(T_all,2));
    T_all(n,:) = T0(I_n); clear I_n
end
clear T0 n

%% 

croptions.pm_labels = "BAC";
croptions.posclass = options.posclass;

tic;
for k = 1:options.KO
    % Load model
    O = struct;
    load(fnames{k},'Omdl');

    O.I_tr = out.O.I_TR{k}; O.X_tr = X(O.I_tr,:);
    O.I_te = out.O.I_TE{k}; O.X_te = X(O.I_te,:);
    
    O.i_keep = out.O.I_KEEP{k}; O.X_te_rd = O.X_te(:,O.i_keep);
    
    for n = 1:N_perm
        O.n.T = T_all(n,:);
        
        O.n.T_tr = logical(O.n.T(O.I_tr));
        O.n.T_te = logical(O.n.T(O.I_te));

        % acquire predictions from RF model
        [O.n.T_pr,O.n.sc] = predict(Omdl,O.X_te_rd);
        O.n.T_pr = str2double(string(cell2mat(O.n.T_pr)));
        O.n.sc = O.n.sc(:,1);

        % compute performance measures
        [O.n.oc,~,~] = ...
            cognemo_classrates(O.n.sc,O.n.T_pr,O.n.T_te,O.I_te,croptions);
        OC_all(n,k) = O.n.oc;
        
        % print progress message
        O.n.it = (k-1)*N_perm + n;
        O.n.percent = round(100*O.n.it/(options.KO*N_perm),2);
        O.n.t2 = toc; O.n.time = datestr(O.n.t2/(24*60*60),"HH:MM:SS");
        O.n.rate = round(O.n.it/O.n.t2,2);
        clc
        fprintf(string(O.n.percent)+"%% done | time: "+...
                string(O.n.time)+" | rate: "+...
                string(O.n.rate)+" iterations/s \n")
        
        O = rmfield(O,'n');
    end
    clear O Omdl
end
clear n k

OCav = mean(OC_all,2);
N_better = length(find(OCav > out.O.OCav));
p = (N_better + 1) / (N_perm + 1);

%%
cd ../..

%% Filter Reliability

% From each outer model, find how many folds each feature occurred in
% I_KEEP

EC.class.rel = zeros([1,EC.in.N_c]);
for k = 1:options.KO
    EC.class.rel(out.O.I_KEEP{1,k}) = EC.class.rel(out.O.I_KEEP{1,k}) + 1;
end
EC.class.rel = EC.class.rel/options.KO;

%% Feature Importance

%% load outer fold models
cd("NCVRF_18-Dec-2021 16-10-29/outermdl")

%%
fnames = strsplit(ls('*.mat'));
fnames = sort(fnames);
fnames = fnames(2:length(fnames));

%% Feature Importance

EC.class.IMP = zeros([options.KO,EC.in.N_c]);
for k = 1:options.KO
    % Load model
    load(fnames{k},'Omdl');
    EC.class.IMP(k,out.O.I_KEEP{1,k}) = Omdl.OOBPermutedPredictorDeltaError;
    clear Omdl
end
cd ../..

%% Rank features by importance

EC.class.MXI = mean(EC.class.IMP,1); EC.class.SXI = std(EC.class.IMP,1);
EC.class.MXI_sc = EC.class.MXI; EC.class.MXI_sc(EC.in.ind_nsc) = 0;
EC.class.MXI_nsc = EC.class.MXI; EC.class.MXI_nsc(EC.in.ind_sc) = 0;
[~,EC.class.top_ind_MXI] = sort(abs(EC.class.MXI),'descend');
[~,EC.class.top_ind_MXI_sc] = sort(abs(EC.class.MXI_sc),'descend');
[~,EC.class.top_ind_MXI_nsc] = sort(abs(EC.class.MXI_nsc),'descend');

%% cgp feature importance (MXI) prep
cg.feat.data = zeros(size(EC.class.MXI));
cg.feat.data(EC.class.top_ind_MXI_nsc(1:pset.cgset.N_top)) = ...
    EC.class.MXI_nsc(EC.class.top_ind_MXI_nsc(1:pset.cgset.N_top));
cg.feat.data = cognemo_symmtx(cg.feat.data,2);
cg.feat.data = cognemo_cg_prep(cg.feat.data,pset.cgset.rlabel.ind_sort);

%% cgp feature importance (MXI) output plot
cg.feat.fig = figure('Units',pset.units,...
                     'Visible',pset.disp,'Position',pset.cgset.figpos);
circularGraph(cg.feat.data,'Label',pset.cgset.labelblank,'Areacolors',pset.cgset.colormap);
cognemo_ringlabel(pset.cgset.rlabel.fnlabsh_sort,pset.cgset.ringoptions)
if pset.save; cg.feat.fname = 'cg_feat.png'; cd ..; cd Figures
saveas(cg.feat.fig,cg.feat.fname); cd ..; cd Pipeline; end

%% bfp feat (MXI) prep
bf.feat.data.ind = EC.class.top_ind_MXI(1:pset.bfset.N_top);
bf.feat.data.MXI = ...
    EC.class.MXI(bf.feat.data.ind);
bf.feat.data.SXI = ...
    EC.class.SXI(bf.feat.data.ind);
%% bfp feat (MXI) plot
bf.feat.fig = figure('Units',pset.units,...
                      'Visible',pset.disp,'Position',pset.bfset.figpos);
cognemo_visualize_feats(bf.feat.data.MXI,...
                        bf.feat.data.SXI,...
                        bf.feat.data.ind,...
                        pset.bfset.options); clear options
set(gca,'FontSize',pset.bfset.fontsize,'TickDir','in')
title(pset.featlabel,'Color','w','FontSize',1.25*get(gca,'FontSize'))
xlabel('Feature Importance')
if pset.save; bf.feat.fname = 'bf_feat_MDD-HC.png'; cd ..; cd Figures
saveas(bf.feat.fig,bf.feat.fname); cd ..; cd Pipeline; end
%% bfp featsc (MXI_sc) prep
bf.featsc.data.ind = EC.class.top_ind_MXI_sc(1:pset.bfset.N_top);
bf.featsc.data.MXI = ...
    EC.class.MXI(bf.featsc.data.ind);
bf.featsc.data.SXI = ...
    EC.class.SXI(bf.featsc.data.ind);
%% bfp featsc (MXI_sc) plot
bf.featsc.fig = figure('Units',pset.units,...
                        'Visible',pset.disp,'Position',pset.bfset.figpos);
cognemo_visualize_feats(bf.featsc.data.MXI,...
                        bf.featsc.data.SXI,...
                        bf.featsc.data.ind,...
                        pset.bfset.options); clear options
set(gca,'FontSize',pset.bfset.fontsize,'TickDir','in')
title(pset.featsclabel,'Color','w','FontSize',1.25*get(gca,'FontSize'))
xlabel('MDD < HC \leftarrow          \rightarrow MDD > HC')
if pset.save; bf.featsc.fname = 'bf_featsc_MDD-HC.png'; cd ..; cd Figures
saveas(bf.featsc.fig,bf.featsc.fname); cd ..; cd Pipeline; end

%% bmp feat (MXI) prep
bm.feat.data.ind = EC.class.top_ind_MXI(1:pset.bmset.N_top);
bm.feat.data.MXI = ...
    EC.class.MXI(bm.feat.data.ind);

%% bmp feat (MXI) plot
bm.feat.fig = figure('Units',pset.units,...
                      'Visible',pset.disp,'Position',pset.bmset.figpos,...
                      'inverthardcopy','off','color','k');
cognemo_areabar(bm.feat.data.MXI,bm.feat.data.ind,pset.bmset.options);
xlabel('Feature Importance')
set(gca,'Position',[0.075 0.125 0.5375 0.75],...
        'FontSize',pset.bmset.fontsize);
title(pset.featlabel,'Color','w','FontSize',1.25*get(gca,'FontSize'))
% annotation('textbox',[0.0125 0.025 0.6375 0.95],...
%            'String',"",...
%            'Color','w','EdgeColor','w','FontSize',24,'LineWidth',2);
cognemo_areabar_legend([0.65 0.15 0.4 0.75],pset.bmset.legoptions)
if pset.save; bm.feat.fname = 'bm_feat_MDD-HC.png'; cd ..; cd Figures
saveas(bm.feat.fig,bm.feat.fname); cd ..; cd Pipeline; end

%% bmp featsc (MXI_sc) prep
bm.featsc.data.ind = EC.class.top_ind_MXI_sc(1:pset.bmset.N_top);
bm.featsc.data.MXI = ...
    EC.class.MXI(bm.featsc.data.ind);

%% bmp featsc (MXI_sc) plot
bm.featsc.fig = figure('Units',pset.units,...
                        'Visible',pset.disp,'Position',pset.bmset.figpos,...
                        'inverthardcopy','off','color','k');
cognemo_areabar(bm.featsc.data.MXI,bm.featsc.data.ind,pset.bmset.options);
xlabel('Feature Importance')
set(gca,'Position',[0.075 0.125 0.5375 0.75],...
        'FontSize',pset.bmset.fontsize);
title(pset.featsclabel,'Color','w','FontSize',1.25*get(gca,'FontSize'))
% annotation('textbox',[0.0125 0.025 0.6375 0.95],...
%            'String',"",...
%            'Color','w','EdgeColor','w','FontSize',24,'LineWidth',2);
cognemo_areabar_legend([0.65 0.15 0.4 0.75],pset.bmset.legoptions)
if pset.save; bm.featsc.fname = 'bm_featsc_MDD-HC.png'; cd ..; cd Figures
saveas(bm.featsc.fig,bm.featsc.fname); cd ..; cd Pipeline; end

%% Network method

%% Load output
cd Output
cd('FNRF_28-Jan-2022 12-17-30')
load('FNRF_28-Jan-2022 12-17-30.mat')
cd ../..

%% image plot

nm.im.data = reshape(out.PMav,[size(out.PMav,1),size(out.PMav,3)]);
nm.im.fig = figure('Units',pset.units,...
                  'Visible',pset.disp,...
                  'OuterPosition',pset.nmset.figpos);
imagesc(nm.im.data)
colormap(pset.nmset.colormap);
colorbar('Color','w','LineWidth',2);
set(gcf,'Color','k','InvertHardCopy','off');
set(gca,'XColor','w','YColor','w',...
        'FontSize',pset.nmset.fontsize,...
        'LineWidth',pset.nmset.linewidth)
xticks(1:size(nm.im.data,2)); xticklabels(out.pm_labels)
yticks(1:size(nm.im.data,1)); yticklabels(out.fnincl');
if pset.save; nm.im.fname = 'networkmethod_allpm.png'; cd ..; cd Figures
saveas(nm.im.fig,nm.im.fname); cd ..; cd Pipeline; end

%% bar plot

nm.ba.x    = 1:size(out.PMav,1);
nm.ba.data = reshape(out.PMav(:,1,4),[1,size(out.PMav,1)]);
nm.ba.std  = reshape(std(out.PM(:,:,4),0,2),[1,size(out.PMav,1)]);
nm.ba.fig = figure('Units',pset.units,...
                   'Visible',pset.disp,...
                   'OuterPosition',pset.nmset.figpos,...
                   'inverthardcopy','off','Color','k');
hold on
nm.ba.bar = barh(nm.ba.x,flip(nm.ba.data), ...
              'FaceColor',pset.pcolor,...
              'EdgeColor','k',...
              'LineStyle','none');
nm.ba.err = errorbar(flip(nm.ba.data),nm.ba.x,flip(nm.ba.std),...
                     'horizontal',...
                     'Marker','.',...
                     'LineWidth',1,...
                     'LineStyle','none',...
                     'Color','w');
set(gca,'color','k','box','on','LineWidth',1,...
        'ytick',[],'xcolor','w','ycolor','w','TickDir','in',...
        'FontName','Helvetica','XColor','w','YColor','w');
yticks(1:size(nm.ba.data,2)); yticklabels(flip(out.fnincl)); xlabel('BAC');
if pset.save; nm.ba.fname = 'networkmethod_bac.png'; cd ..; cd Figures
saveas(nm.ba.fig,nm.ba.fname); cd ..; cd Pipeline; end
%}
%% BAC for different random seeds, for different networks

load('FNRF_grand.mat')
cd Output/FNRF_Feb23
fnames = strsplit(ls); cd([fnames{1} ' ' fnames{2}]);
fname = strtrim(ls('*.mat')); load(fname);
fnlabels = out.fnincl; clear options out
pm.fig = figure('Units',pset.units,...
                'Visible','on',...
                'OuterPosition',[0 0 15 8],...
                'inverthardcopy','off','Color','k');
X = 1:size(PM_grand,4);
hold on
for i = 1:size(PM_grand,1)
    Yi = reshape(mean(PM_grand(i,4,:,:),3),[1,size(PM_grand,4)]);
    Ei = reshape(std(PM_grand(i,4,:,:),[],3),[1,size(PM_grand,4)]);
    err_i = errorbar(Yi,Ei,'Color',pset.pcolor,'LineStyle','none');
end
for i = 1:size(PM_grand,1)
    Yi = reshape(mean(PM_grand(i,4,:,:),3),[1,size(PM_grand,4)]);
    line_i = plot(X,Yi,'Color','r','Marker','.');
end
Y = mean(PM_grand(:,4,:,:),3);
Y = reshape(mean(Y,1),[1,size(PM_grand,4)]);
line = plot(X,Y,'MarkerFaceColor','b','Color','w','Marker','o','LineWidth',2);
Ymax = 0:(size(PM_grand,4)+1);
linemax = plot(Ymax,max(Y)*ones(size(Ymax)),'Marker','none','Color','g','LineWidth',2);
set(gca,'color','k','box','on','LineWidth',1,...
        'xcolor','w','ycolor','w','TickDir','in',...
        'FontName','Helvetica','XColor','w','YColor','w','fontsize',18);
xlim([0,(size(PM_grand,4)+1)])
xlabels = unique(pset.cgset.rlabel.fns.og,'stable');
xticks(1:size(PM_grand,4))
xticklabels(xlabels); clear xlabels
ylabels = 0.1*(0:9); ylabels = [ylabels round(max(Y),4)]; ylabels = sort(ylabels,'ascend');
yticks(ylabels)
yticklabels(string(ylabels))
cd ../../..

%% Feature Importance---one run
%{
% DMN index = 6;
dmn.i = 6;
dmn.I = logical(out.fs.I(dmn.i,:));

dmn.MXI = out.imp.mean(dmn.I,1,dmn.i);
dmn.SXI = out.imp.std(dmn.I,1,dmn.i);
[~,dmn.I_sort]  = sort(dmn.MXI,'descend');

%% bfp feat (MXI) prep
bf.feat.data.ind = dmn.I_sort(1:pset.bfset.N_top);
bf.feat.data.MXI = dmn.MXI(bf.feat.data.ind);
bf.feat.data.SXI = dmn.SXI(bf.feat.data.ind);

%% bfp feat (MXI) plot
bf.feat.fig = figure('Units',pset.units,...
                      'Visible',pset.disp,'Position',pset.bfset.figpos);
cognemo_visualize_feats(bf.feat.data.MXI,...
                        bf.feat.data.SXI,...
                        bf.feat.data.ind,...
                        pset.bfset.options); clear options
set(gca,'FontSize',pset.bfset.fontsize,'TickDir','in')
title(pset.featlabel,'Color','w','FontSize',1.25*get(gca,'FontSize'))
xlabel('Feature Importance')
if pset.save; bf.feat.fname = 'bf_feat_MDD-HC.png'; cd ..; cd Figures
saveas(bf.feat.fig,bf.feat.fname); cd ..; cd Pipeline; end

%% bmp feat (MXI) prep
bm.feat.data.ind = dmn.I_sort(1:pset.bmset.N_top);
bm.feat.data.MXI = ...
    dmn.MXI(bm.feat.data.ind);

%% bmp feat (MXI) plot
bm.feat.fig = figure('Units',pset.units,...
                      'Visible',pset.disp,'Position',pset.bmset.figpos,...
                      'inverthardcopy','off','color','k');
cognemo_areabar(bm.feat.data.MXI,bm.feat.data.ind,pset.bmset.options);
xlabel('Feature Importance')
set(gca,'Position',[0.075 0.125 0.5375 0.75],...
        'FontSize',pset.bmset.fontsize);
title(pset.featlabel,'Color','w','FontSize',1.25*get(gca,'FontSize'))
% annotation('textbox',[0.0125 0.025 0.6375 0.95],...
%            'String',"",...
%            'Color','w','EdgeColor','w','FontSize',24,'LineWidth',2);
cognemo_areabar_legend([0.65 0.15 0.4 0.75],pset.bmset.legoptions)
if pset.save; bm.feat.fname = 'bm_feat_MDD-HC.png'; cd ..; cd Figures
saveas(bm.feat.fig,bm.feat.fname); cd ..; cd Pipeline; end
%}
%% Feature Importance---stability (Horiz)
%{
load('FNRF_grand.mat')

N_top = 15;
N_se  = size(FI_grand,1);

cd Output/FNRF_Nov22
fnames = strsplit(ls); cd([fnames{1} ' ' fnames{2}]);
fname = strtrim(ls('*.mat')); load(fname);
fnlabels = out.fnincl; clear options out
bf.fig = figure('Units',pset.units,...
                'Visible','on',...
                'OuterPosition',[0 0 15 8],...
                'inverthardcopy','off','Color','k');
% Prepare
X = 1:N_top;
FI_ave1 = mean(FI_grand,1); FI_ave2 = mean(FI_ave1,3);
FI_grandave = reshape(FI_ave2,[1,size(FI_ave2,2)]);
[Y,I_top] = maxk(abs(FI_grandave),N_top); clear FI_ave1 FI_ave2

hold on
%{
for i = 1:size(PM_grand,1)
    Yi = reshape(mean(PM_grand(i,4,:,:),3),[1,size(PM_grand,4)]);
    Ei = reshape(std(PM_grand(i,4,:,:),[],3),[1,size(PM_grand,4)]);
    err_i = errorbar(Yi,Ei,'Color',pset.pcolor,'LineStyle','none');
end
%}
for i = 1:N_se
    FI_ave1_i = reshape(mean(FI_grand(i,:,:),3),[1,size(FI_grand,2)]);
    Yi = FI_ave1_i(I_top);
    line_i = plot(X,Yi,'Color','r','Marker','.');
end

line = plot(X,Y,'MarkerFaceColor','b','Color','w','Marker','o','LineWidth',2);
Xmax = 0:(N_top+1);
linemax = plot(Xmax,max(Y)*ones(size(Xmax)),'Marker','none','Color','g','LineWidth',2);
set(gca,'color','k','box','on','LineWidth',1,...
        'xcolor','w','ycolor','w','TickDir','in',...
        'FontName','Helvetica','XColor','w','YColor','w','fontsize',18);
xlim([0,(N_top+1)])

xlabels = pset.bfset.options.clabeldir.og(I_top);

xticks(1:N_top)
xticklabels(xlabels); clear xlabels
ylabels = 0.1*(0:9); ylabels = [ylabels max(Y)]; ylabels = sort(ylabels,'ascend');
yticks(ylabels)
yticklabels(string(ylabels))
cd ../../..
%}
%% Feature Importance---stability (Vert)

N_top = 20;

dmn.i = 6;

load('FNRF_grand.mat')
cd Output/FNRF_Feb23

fnames = strsplit(ls); cd([fnames{1} ' ' fnames{2}]);
fname = strtrim(ls('*.mat')); load(fname);

fnlabels = out.fnincl; 
dmn.I = find(out.fs.I(dmn.i,:));
N_se  = size(FI_grand,1);

clear options out

outerpos = [0 0 18 12];
bf.fig = figure('Units',pset.units,...
                'Visible','on',...
                'OuterPosition',[0 0 15 8],...
                'inverthardcopy','off','Color','k');
% Prepare
X = 1:N_top;
FI_ave1 = mean(FI_grand,1); FI_ave2 = mean(FI_ave1,3);
FI_grandave = reshape(FI_ave2,[1,size(FI_ave2,2)]);
[Y,I_top] = maxk(abs(FI_grandave),N_top); clear FI_ave1 FI_ave2
Y = flip(Y);

hold on
%{
for i = 1:size(PM_grand,1)
    Yi = reshape(mean(PM_grand(i,4,:,:),3),[1,size(PM_grand,4)]);
    Ei = reshape(std(PM_grand(i,4,:,:),[],3),[1,size(PM_grand,4)]);
    err_i = errorbar(Yi,Ei,'Color',pset.pcolor,'LineStyle','none');
end
%}
for i = 1:N_se
    FI_ave1_i = reshape(mean(FI_grand(i,:,:),3),[1,size(FI_grand,2)]);
    Yi = FI_ave1_i(I_top);
    
    Yi = flip(Yi);
    
    line_i = plot(Yi,X,'Color','r','Marker','.');
end

line = plot(Y,X,'MarkerFaceColor','b','Color','w','Marker','o','LineWidth',2);
Ymax = 0:(N_top+1);
%linemax = plot(Ymax,max(X)*ones(size(Ymax)),'Marker','none','Color','g','LineWidth',2);
set(gca,'color','k','box','on','LineWidth',1,...
        'xcolor','w','ycolor','w','TickDir','in',...
        'FontName','Helvetica','XColor','w','YColor','w','fontsize',14*outerpos(4)/12);
ylim([0,(N_top+1)])

ylabels = flip(pset.bfset.options.clabeldir.og(dmn.I(I_top)));

yticks(1:N_top)
yticklabels(ylabels); clear ylabels outerpos
xlabel('Feature Importance')

cd ../../..
%}

%% MADRS correlations: Data

load("data.mat")
load('FNRF_grand.mat')
cd Output/FNRF_Feb23

fnames = strsplit(ls); cd([fnames{1} ' ' fnames{2}]);
fname = strtrim(ls('*.mat')); load(fname); clear fname fnames; cd ../../..

dmn.i   = 6;
dmn.I   = find(out.fs.I(dmn.i,:));
dmn.N_c = length(dmn.I);
dmn.labels = pset.bfset.options.clabeldir.og(dmn.I);

mcor.Y = EC.in.X(:,dmn.I);  % Effective connectivity values
mcor.X = data.V(:,1);          % MADRS scores

mcor.rho = zeros(1,dmn.N_c);
mcor.p   = zeros(1,dmn.N_c);

for j = 1:dmn.N_c
    [mcor.rho(j),mcor.p(j)] = corr(mcor.X,mcor.Y(:,j));
end
clear j

mcor.J = find(mcor.p<0.05); mcor.N_J = length(mcor.J);

mcor.msg = "The MADRS sum scores correlate with " + ...
              string(mcor.N_J) + ...
              " within-DMN connections\n";
fprintf(mcor.msg)

rho  = mcor.rho(mcor.J)';
p    = mcor.p(mcor.J)';
label = dmn.labels(mcor.J)';
mcor_table = table(label,rho,p); clear label rho p
mcor_name  = 'mcor_table.csv';
writetable(mcor_table,mcor_name)

%% T-stat analysis---DMN
dmn.i = 6;
dmn.I = find(out.fs.I(dmn.i,:));
N_top = 15; foptions.N_top = N_top;

load('in.mat')
EC.rdata.X = X; data.T = T;
xset.toptions.corr = "none";

EC.pdata = cognemo_prepC(EC,data);
EC = cognemo_mass(EC,xset);

DMN.TXD  = -EC.tstat.TXD_f(dmn.I); DMN.SXD = EC.tstat.SXD_f(dmn.I);
indp = (DMN.TXD > 0); indn = (DMN.TXD < 0);
DMN.TXDp = DMN.TXD(indp); DMN.SXDp = DMN.SXD(indp); clear indp
DMN.TXDn = DMN.TXD(indn); DMN.SXDn = DMN.SXD(indn); clear indn

[~,DMN.top_I]  = maxk(abs(DMN.TXD),N_top);
[~,DMN.top_Ip] = maxk(abs(DMN.TXDp),N_top);
[~,DMN.top_In] = maxk(abs(DMN.TXDn),N_top);

DMN.SXDtop = DMN.SXD(DMN.top_I);
DMN.TXDtop = DMN.TXD(DMN.top_I);
DMN.I_top  = dmn.I(DMN.top_I);

DMN.SXDptop = DMN.SXDp(DMN.top_Ip);
DMN.TXDptop = DMN.TXDp(DMN.top_Ip);
DMN.Ip_top  = dmn.I(DMN.top_Ip);

DMN.SXDntop = DMN.SXDn(DMN.top_In);
DMN.TXDntop = DMN.TXDn(DMN.top_In);
DMN.In_top  = dmn.I(DMN.top_In);

xset.foptions.clabel    = pset.bfset.options.clabel;
xset.foptions.dir       = 1;
xset.foptions.clabeldir = pset.bfset.options.clabeldir;
xset.foptions.xlabel    = "t-value of MDD-HC";

figure('Units',pset.units,...
       'Visible','on',...
       'OuterPosition',[0 0 7 8],...
       'inverthardcopy','off','Color','k');
cognemo_visualize_feats(DMN.TXDtop,...
                        DMN.SXDtop,...
                        DMN.I_top,...
                        xset.foptions);

figure('Units',pset.units,...
       'Visible','on',...
       'OuterPosition',[0 0 7 8],...
       'inverthardcopy','off','Color','k');
cognemo_visualize_feats(DMN.TXDptop,...
                        DMN.SXDptop,...
                        DMN.Ip_top,...
                        xset.foptions);
                    
figure('Units',pset.units,...
       'Visible','on',...
       'OuterPosition',[0 0 7 8],...
       'inverthardcopy','off','Color','k');
cognemo_visualize_feats(DMN.TXDntop,...
                        DMN.SXDntop,...
                        DMN.In_top,...
                        xset.foptions);

%% T-stat analysis---DMN self-connections
%
N_c = N_r*N_r;
ind_sc = find(reshape(eye(N_r),[1,N_c]));
dmn.I_sc = intersect(dmn.I,ind_sc); clear ind_sc

DMN.TXD_sc  = EC.tstat.TXD_f(dmn.I_sc);
DMN.SXD_sc = EC.tstat.SXD_f(dmn.I_sc);
[~,DMN.top_I_sc]  = maxk(abs(DMN.TXD_sc),N_top);

DMN.SXDtop_sc = DMN.SXD_sc(DMN.top_I_sc);
DMN.TXDtop_sc = DMN.TXD_sc(DMN.top_I_sc);
DMN.I_sc_top  = dmn.I_sc(DMN.top_I_sc);

xset.foptions.clabel    = pset.bfset.options.clabel;
xset.foptions.dir       = 1;
xset.foptions.clabeldir = pset.bfset.options.clabeldir;

figure('Units',pset.units,...
       'Visible','on',...
       'OuterPosition',[0 0 7 8],...
       'inverthardcopy','off','Color','k');
cognemo_visualize_feats(DMN.TXDtop_sc,...
                        DMN.SXDtop_sc,...
                        DMN.I_sc_top,...
                        xset.foptions);
%}

%% T-stat analysis---DMN non-self-connections

ind_sc = find(reshape(eye(N_r),[1,N_c]));
dmn.I_nsc = setdiff(dmn.I,ind_sc); clear ind_sc

DMN.MXD_nsc = EC.tstat.MXD_f(dmn.I_nsc);
DMN.TXD_nsc = -EC.tstat.TXD_f(dmn.I_nsc);
DMN.SXD_nsc = EC.tstat.SXD_f(dmn.I_nsc);
[~,DMN.top_I_nsc]  = maxk(abs(DMN.TXD_nsc),N_top);

DMN.MXDtop_nsc = DMN.MXD_nsc(DMN.top_I_nsc);
DMN.SXDtop_nsc = DMN.SXD_nsc(DMN.top_I_nsc);
DMN.TXDtop_nsc = DMN.TXD_nsc(DMN.top_I_nsc);
DMN.I_nsc_top  = dmn.I_nsc(DMN.top_I_nsc);

xset.foptions.clabel    = pset.bfset.options.clabel;
xset.foptions.dir       = 1;
xset.foptions.clabeldir = pset.bfset.options.clabeldir;

figure('Units',pset.units,...
       'Visible','on',...
       'OuterPosition',[0 0 7 8],...
       'inverthardcopy','off','Color','k');
cognemo_visualize_feats(DMN.TXDtop_nsc,...
                        DMN.SXDtop_nsc,...
                        DMN.I_nsc_top,...
                        xset.foptions);


%%

DMN.TXD_nsc_ave = mean(DMN.TXD_nsc);
if DMN.TXD_nsc_ave > 0
    DMN.avemsg = ...
        "In MDD, average within-DMN excitatory connectivity INCREASED\n";
elseif DMN.TXD_nsc_ave < 0
    DMN.avemsg = ...
        "In MDD, average within-DMN excitatory connectivity DECREASED\n";
else
    DMN.avemsg = ...
        "In MDD, average within-DMN excitatory connectivity DID NOT CHANGE]n";
end

fprintf(DMN.avemsg)

%% PM value table

PM_grandave = mean(PM_grand,1);
PM_grandave = mean(PM_grandave,3);
PM_grandave = reshape(PM_grandave,[7,15]);
PM_grandave = PM_grandave';

PM_labels = {'ACC','SE','SP','BAC','PPV','NPV','AUC'};

[fn_sort,~,i_fn_sort] = unique(options.fnlabel);
fn_count = accumarray(i_fn_sort,1).^2; clear i_fn_sort
[~,I_fn_sort] = sort(unique(options.fnlabel,'stable'),'ascend');
PM_table = round(PM_grandave,4); PM_table = PM_table(I_fn_sort,:);

tabcat1 = table(fn_sort,'VariableNames',{'Network'});
tabcat2 = table(fn_count,'VariableNames',{'Connections'});
tabcat3 = array2table(PM_table,'VariableNames',PM_labels);

PM_table = [ tabcat1 tabcat2 tabcat3 ];

PMtab_name  = 'PM_table.csv';
writetable(PM_table,PMtab_name)

%% ROC plots

figure('Units',pset.units,...
       'Visible','on',...
       'OuterPosition',[0 0 15 9],...
       'inverthardcopy','off','Color','k');
N_sc = 5; N_sr = 3;
for k = 1:N_sc
    for j = 1:N_sr
        i_ax = (j-1)*N_sc + k; i_roc = I_fn_sort(i_ax);
        subaxis(N_sr,N_sc,i_ax,'SpacingHoriz',0.01)
        hold on
        for i = 1:50
            plot(ROCx_grand{i,i_roc},ROCy_grand{i,i_roc},'Color','r','LineWidth',0.5)
        end
        % title(out.fnincl(i_roc),'FontSize',10,'Color','w')
        text(0.98,1.02,out.fnincl(i_roc),...
             'HorizontalAlignment','Right','VerticalAlignment','Bottom',...
             'FontSize',14,'Color','w','FontWeight','Bold');
        set(gca,'color','k','box','on','LineWidth',1,...
        'xcolor','w','ycolor','w','TickDir','in','DataAspectRatio',[1,1,1],...
        'FontName','Helvetica','XColor','w','YColor','w','fontsize',14);
        % auc_dim = [0.5 0.5 0.1 0.1];
        auc_str = ['$\overline{AUC}=' char(string(round(PM_grandave(i_roc,7),4))) '$'];
        text(0.25,0.12,auc_str,'Color','w','FontSize',14,'Interpreter','latex');
        if k ~= 1
            set(gca,'YTick',[]);
        else
            ylabel('True positive rate','FontSize',12);
        end
        if j ~= N_sr
            set(gca,'XTick',[]);
        else
            xlabel('False positive rate','FontSize',12);
        end
        hold off
        
    end
end

%% Feature collinearity

% covariate correction
load('data.mat'); % load covariate data
EC.dmn = dmn;
EC.dmn.X = EC.pdata.X(:,EC.dmn.I);

[~,EC.dmn.E] = cognemo_cctrain(EC.dmn.X,data.V,"none");

% correlations
CC = corrcoef(EC.dmn.E); CC(logical(eye(length(CC)))) = 0;
CC_vec = abs(reshape(CC,[1,length(CC)^2]));
hist.data = CC_vec(logical(CC_vec)); clear CC CC_vec

%% Correlation Coefficients Graph - EC plot

hist.set.barcolor = 'm';
hist.vline.set.x1 = 0.3;
hist.vline.set.color1 = 'g';

hist.vline.set.x2 = 0.5;
hist.vline.set.color2 = 'b';

hist.vline.data.y1 = 100 * (length(find(hist.data >= hist.vline.set.x1)) / ...
                            length(hist.data));
hist.vline.data.y2 = 100 * (length(find(hist.data >= hist.vline.set.x2)) / ...
                            length(hist.data));
                        
hist.vline.set.label1 = string(round(hist.vline.data.y1,2,'significant')) + ...
                        "$\%$ of correlations $|r| \geq$ " + string(hist.vline.set.x1);
hist.vline.set.label2 = string(round(hist.vline.data.y2,2,'significant')) + ...
                        "$\%$ of correlations $|r| \geq$ " + string(hist.vline.set.x2);

% Histogram
figure('OuterPosition',[50 50 600 600],...
       'InvertHardCopy','off','Color','k')
                    
hist.draw = histogram(hist.data,'DisplayStyle','stairs',...
                                'EdgeColor',hist.set.barcolor,...
                                'LineWidth',1.5);
                            
hist.vline.draw1 = line([hist.vline.set.x1,hist.vline.set.x1],...
                        get(gca,'Ylim'),...
                        'Color',hist.vline.set.color1,'LineWidth',1);
hist.vline.draw2 = line([hist.vline.set.x2,hist.vline.set.x2],...
                        get(gca,'Ylim'),...
                        'Color',hist.vline.set.color2,'LineWidth',1);
                            
set(gca,'box','off','Color','k','XColor','w','YColor','w',...
        'LineWidth',1.5,'fontsize',14,...
        'Xgrid','on','Ygrid','on','GridLineStyle',':','GridAlpha',0.2);                
    
xticks(0.1*(0:10)); xlim([0,1]);
xlabel('Absolute correlation coefficient')
ylabel('Number of feature correlations')

legend([hist.vline.draw1,hist.vline.draw2],...
       [hist.vline.set.label1,hist.vline.set.label2],...
       'FontSize',18,'TextColor','w','Interpreter','latex');

%% DMN-only connectograms
%{

Requirements:
- Select DMN connection indices
- Select data and labels
- rearrange indices for cgp if needed

%}

% labels
if isfield(dmn,"I") dmn.I_c = dmn.I; dmn = rmfield(dmn,"I"); end
if isfield(dmn,"labels") dmn.clabels = dmn.labels; dmn = rmfield(dmn,"labels"); end

dmn.I_r = find(pset.cgset.rlabel.fns.og=="DMN");
dmn.rlabels = pset.cgset.rlabel.rlabel.og(dmn.I_r);
EC.dmn = dmn;
EC.dmn.MX0 = zeros(size(EC.stats.MX0));
EC.dmn.MX1 = zeros(size(EC.stats.MX1));
EC.dmn.TXD = zeros(size(EC.stats.TXD));

EC.dmn.MX0(EC.dmn.I_c) = EC.stats.MX0(EC.dmn.I_c);
EC.dmn.MX1(EC.dmn.I_c) = EC.stats.MX1(EC.dmn.I_c);
EC.dmn.TXD(EC.dmn.I_c) = EC.stats.TXD(EC.dmn.I_c);


%% cgp meanHC (MX0) prep
EC.dmn.MX0_sc  = EC.dmn.MX0; EC.dmn.MX0_sc(EC.in.ind_nsc) = 0;
EC.dmn.MX0_nsc = EC.dmn.MX0; EC.dmn.MX0_nsc(EC.in.ind_sc) = 0;
[~,EC.dmn.top_ind_MX0] = sort(abs(EC.dmn.MX0),'descend');
[~,EC.dmn.top_ind_MX0_sc] = sort(abs(EC.dmn.MX0_sc),'descend');
[~,EC.dmn.top_ind_MX0_nsc] = sort(abs(EC.dmn.MX0_nsc),'descend');

datatemp = zeros(size(EC.dmn.MX0));
datatemp(EC.dmn.top_ind_MX0_nsc(1:pset.cgset.N_top)) = ...
    EC.dmn.MX0(EC.dmn.top_ind_MX0_nsc(1:pset.cgset.N_top));
cg.dmn.meanHC.data = datatemp(EC.dmn.I_c);
cg.dmn.meanHC.data = cognemo_symmtx(cg.dmn.meanHC.data,2);
cg.dmn.meanHC.data = cognemo_cg_prep(cg.dmn.meanHC.data,1:numel(EC.dmn.I_r));
% change area colors
colormaptemp = pset.cgset.fn.colormap(pset.cgset.rlabel.fn_unsort,:);
cg.dmn.colormap = colormaptemp(EC.dmn.I_r,:);
% labels
cg.dmn.labels = dmn.rlabels;
for i = 1:numel(dmn.rlabels)
    cg.dmn.rlabels(i) = erase(dmn.rlabels(i),"_");
end
cg.dmn.labels = cellstr(cg.dmn.rlabels);



%% cgp meanHC (MX0) output plot
cg.dmn.meanHC.fig = figure('Units',pset.units,...
                       'Visible',pset.disp,'Position',pset.cgset.figpos);
circularGraph(cg.dmn.meanHC.data,'Label',cg.dmn.labels,'Areacolors',cg.dmn.colormap);
% cognemo_ringlabel(pset.cgset.rlabel.fns.og(pset.cgset.rlabel.fn_sort),pset.cgset.fn.ringoptions)
%{
if pset.save; cg.dmn.meanHC.fname = 'cgfn_meanHC.png'; cd ..; cd Figures
saveas(cg.meanHC.fig,cg.meanHC.fname); cd ..; cd Pipeline; end
%}

%% cgp meanHC (MX1) prep
EC.dmn.MX0_sc  = EC.dmn.MX1; EC.dmn.MX1_sc(EC.in.ind_nsc) = 0;
EC.dmn.MX1_nsc = EC.dmn.MX1; EC.dmn.MX1_nsc(EC.in.ind_sc) = 0;
[~,EC.dmn.top_ind_MX1] = sort(abs(EC.dmn.MX1),'descend');
[~,EC.dmn.top_ind_MX1_sc] = sort(abs(EC.dmn.MX1_sc),'descend');
[~,EC.dmn.top_ind_MX1_nsc] = sort(abs(EC.dmn.MX1_nsc),'descend');

datatemp = zeros(size(EC.dmn.MX1));
datatemp(EC.dmn.top_ind_MX1_nsc(1:pset.cgset.N_top)) = ...
    EC.dmn.MX1(EC.dmn.top_ind_MX1_nsc(1:pset.cgset.N_top));
cg.dmn.meanHC.data = datatemp(EC.dmn.I_c);
cg.dmn.meanHC.data = cognemo_symmtx(cg.dmn.meanHC.data,2);
cg.dmn.meanHC.data = cognemo_cg_prep(cg.dmn.meanHC.data,1:numel(EC.dmn.I_r));
% change area colors
colormaptemp = pset.cgset.fn.colormap(pset.cgset.rlabel.fn_unsort,:);
cg.dmn.colormap = colormaptemp(EC.dmn.I_r,:);
% labels
cg.dmn.labels = dmn.rlabels;
for i = 1:numel(dmn.rlabels)
    cg.dmn.rlabels(i) = erase(dmn.rlabels(i),"_");
end
cg.dmn.labels = cellstr(cg.dmn.rlabels);



%% cgp meanHC (MX1) output plot
cg.dmn.meanHC.fig = figure('Units',pset.units,...
                       'Visible',pset.disp,'Position',pset.cgset.figpos);
circularGraph(cg.dmn.meanHC.data,'Label',cg.dmn.labels,'Areacolors',cg.dmn.colormap);
% cognemo_ringlabel(pset.cgset.rlabel.fns.og(pset.cgset.rlabel.fn_sort),pset.cgset.fn.ringoptions)
%{
if pset.save; cg.dmn.meanHC.fname = 'cgfn_meanHC.png'; cd ..; cd Figures
saveas(cg.meanHC.fig,cg.meanHC.fname); cd ..; cd Pipeline; end
%}

%% cgp meanHC (TXD) prep
EC.dmn.MX0_sc  = EC.dmn.TXD; EC.dmn.TXD_sc(EC.in.ind_nsc) = 0;
EC.dmn.TXD_nsc = EC.dmn.TXD; EC.dmn.TXD_nsc(EC.in.ind_sc) = 0;
[~,EC.dmn.top_ind_TXD] = sort(abs(EC.dmn.TXD),'descend');
[~,EC.dmn.top_ind_TXD_sc] = sort(abs(EC.dmn.TXD_sc),'descend');
[~,EC.dmn.top_ind_TXD_nsc] = sort(abs(EC.dmn.TXD_nsc),'descend');

datatemp = zeros(size(EC.dmn.TXD));
datatemp(EC.dmn.top_ind_TXD_nsc(1:pset.cgset.N_top)) = ...
    EC.dmn.TXD(EC.dmn.top_ind_TXD_nsc(1:pset.cgset.N_top));
cg.dmn.meanHC.data = datatemp(EC.dmn.I_c);
cg.dmn.meanHC.data = cognemo_symmtx(cg.dmn.meanHC.data,2);
cg.dmn.meanHC.data = cognemo_cg_prep(cg.dmn.meanHC.data,1:numel(EC.dmn.I_r));
% change area colors
colormaptemp = pset.cgset.fn.colormap(pset.cgset.rlabel.fn_unsort,:);
cg.dmn.colormap = colormaptemp(EC.dmn.I_r,:);

%% Labels and reordering

DMNtab = readtable('RegionLabels_Apr22.xlsx','Sheet','Sheet1','Range','Q1:R271');
dmnorder = table2array(DMNtab(:,1));
dmnlabel = table2array(DMNtab(:,2));
EC.dmn.cg_ind = dmnorder(EC.dmn.I_r); clear dmnorder DMNtab
labeltemp = dmnlabel(EC.dmn.I_r); clear dmnlabel
for i = 1:numel(labeltemp)
    % cg.dmn.rlabels(i) = i;
    labeltemp(i) = erase(labeltemp(i),"_");
end
[~,EC.dmn.cg_ind] = sort(EC.dmn.cg_ind,'ascend');


EC.dmn.N_r = sqrt(EC.dmn.N_c);
cg_rotate_ind = 1:EC.dmn.N_r; ind_temp = cg_rotate_ind;
cg_rotate_ind(1:0.25*EC.dmn.N_r) = ind_temp(0.75*EC.dmn.N_r+1:EC.dmn.N_r);
cg_rotate_ind(0.25*EC.dmn.N_r+1:EC.dmn.N_r) = ind_temp(1:0.75*EC.dmn.N_r);
EC.dmn.cg_ind = EC.dmn.cg_ind(cg_rotate_ind);

cg.dmn.rlabel = labeltemp(EC.dmn.cg_ind); clear labeltemp
cg.dmn.meanHC.data = cg.dmn.meanHC.data(EC.dmn.cg_ind,EC.dmn.cg_ind);

cg.dmn.labels = cellstr(cg.dmn.rlabel);

%% cgp meanHC (TXD) output plot
cg.dmn.meanHC.fig = figure('Units',pset.units,...
                           'Visible',pset.disp,'Position',[0 0 12 12]);
circularGraph(cg.dmn.meanHC.data,'Label',cg.dmn.labels,'Areacolors',cg.dmn.colormap);

cognemo_alphabar(cg.dmn.meanHC.data,[.1 .12 .06 .25],fig1.set.cbar);


%% Plot mcor in connectogram

EC.dmn.mcor = zeros(size(mcor.rho));
EC.dmn.mcor(mcor.J) = mcor.rho(mcor.J);
[~,ind_top_mcor] = maxk(abs(EC.dmn.mcor),238);

cg.dmn.mcor.data = zeros(size(EC.dmn.mcor));
cg.dmn.mcor.data(ind_top_mcor) = EC.dmn.mcor(ind_top_mcor); clear ind_top_mcor

cg.dmn.mcor.data = cognemo_symmtx(cg.dmn.mcor.data,3);
cg.dmn.mcor.data = cognemo_cg_prep(cg.dmn.mcor.data,1:numel(EC.dmn.I_r));

cg.dmn.mcor.data = cg.dmn.mcor.data(EC.dmn.cg_ind,EC.dmn.cg_ind);

cg.dmn.mcor.fig = figure('Units',pset.units,...
                         'Visible',pset.disp,'Position',[0 0 12 12]);
circularGraph(cg.dmn.mcor.data,'Label',cg.dmn.labels,'Areacolors',cg.dmn.colormap);

cognemo_alphabar(cg.dmn.mcor.data,[.1 .12 .07 .25],fig1.set.cbar);


