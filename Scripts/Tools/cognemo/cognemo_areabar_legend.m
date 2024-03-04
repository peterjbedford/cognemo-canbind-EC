function legax = cognemo_areabar_legend(dim,legoptions)
%% Preamble
%{
%}
%% Unpack options

linespacing = legoptions.linespacing;
fontsize = legoptions.fontsize;
arealabels = legoptions.arealabels;
areacolor = legoptions.areacolor./255;

%% Define Colourmap

areastr = unique(arealabels,'stable');
areacolor = unique(areacolor,'rows','stable');

N_area = length(areastr);

%% Plot legend

figdim = get(gcf,'Position');

ywidth_inches = figdim(4)*dim(4);
linewidth_points = linespacing*fontsize; inches_points = 1/72;
linewidth_inches = linewidth_points*inches_points;
N_lines = floor(ywidth_inches/linewidth_inches);

ylim_u = dim(2)+dim(4); % normalized to figure


%areacolor = flip(areacolor); areastr = flip(areastr);

for j = 1:N_area
    ypos =  ylim_u - dim(4)*(j/N_lines);
    dim(2) = ypos;
    color = char(string(areacolor(j,1)) + "," + ...
                 string(areacolor(j,2)) + "," + ...
                 string(areacolor(j,3)));
    txtcol = ['{\color[rgb]{' color '} \bullet} '];
    txtlab = char(areastr(j));
    txt = [txtcol txtlab];
    annotation('textbox',dim,'String',txt,...
               'LineStyle','none',...
               'VerticalAlignment','bottom',...
               'Color','w','FontSize',fontsize);
end
