function cognemo_visualize_feats(top_mean,top_std,top_ind,foptions)
%% Preamble
%{
Displays performance measure charts for different classification models
%}
%% unpack parameter options

N_top = min([length(top_mean),length(top_std),length(top_ind)]);
if isfield(foptions,'N_top')
    N_top = foptions.N_top;
end

shared_ind = [];
if isfield(foptions,'shared_ind')   
    shared_ind = foptions.shared_ind;
end

barcolor_s = [0 0 0];
if isfield(foptions,'shared_color')
    barcolor_s = foptions.shared_color;
end

clabel = foptions.clabel;
if foptions.dir
    clabel = foptions.clabeldir;
end

%% Prepare data

topN_mean = top_mean(1:N_top);
topN_std  = top_std(1:N_top);
topN_ind  = top_ind(1:N_top);
[~,~,s_ind_plot] = intersect(shared_ind,topN_ind,'stable');
ns_ind_plot = setdiff(1:N_top,s_ind_plot);

%% Plot details

clabel_plot = clabel.og(topN_ind);

%% plot 1: top important features (mean +/- std)

fig1 = gcf; hold on;

y = flip(topN_mean); 
x_s = N_top+1-s_ind_plot;
y_s = zeros(size(y)); y_s(x_s) = y(x_s);
x_ns = N_top+1-ns_ind_plot;
y_ns = zeros(size(y)); y_ns(x_ns) = y(x_ns);

ns_alpha = 0.5;
if isempty(shared_ind)
    ns_alpha = 1;
end

pcolor = [255 181 0]./255; ncolor = [0 201 239]./255;
ind_p  = y_ns > 0; ind_n = y_ns < 0;
y_ns_p = zeros(size(y_ns)); y_ns_p(ind_p) = y_ns(ind_p);
y_ns_n = zeros(size(y_ns)); y_ns_n(ind_n) = y_ns(ind_n);
b_ns_p = barh(1:N_top,y_ns_p,...
              'FaceColor',pcolor,...
              'EdgeColor','k',...
              'LineStyle','none',...
              'FaceAlpha',ns_alpha);
b_ns_n = barh(1:N_top,y_ns_n,...
              'FaceColor',ncolor,...
              'EdgeColor','k',...
              'LineStyle','none',...
              'FaceAlpha',ns_alpha);


if isfield(foptions,'errbar')
    if foptions.errbar == "std" || foptions.errbar == "ste"
        err = flip(topN_std);
        if foptions.errbar == "ste"
            k = foptions.k;
            err = err./sqrt(k);
        end
    end
    e1 = errorbar(y,1:N_top,err,...
            'horizontal',...
            'Marker','.',...
            'LineWidth',1,...
            'LineStyle','none',...
            'Color','w');
end

set(gca,'color','k','box','on','LineWidth',1,...
        'ytick',[],'xcolor','w','ycolor','w','TickDir','in',...
        'FontName','Helvetica','XColor','w','YColor','w');
set(gcf,'inverthardcopy','off','color','k'); 

ax1 = gca;
xlabel(foptions.xlabel)
yrule1 = ax1.YAxis; xrule1 = ax1.XAxis;
yticks(1:N_top)
ax1.YTickLabel = flip(clabel_plot); 

end