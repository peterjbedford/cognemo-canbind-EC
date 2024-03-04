function ax = cognemo_alphabar(data,position,options)
%% Preamble
%{
%}
%% Vertical alphabar

fontsize = options.fontsize;
colormin = options.colormin;
colormax = options.colormax;

% Define max and min
if ~isfield(options,'lims')
    ymin = round(min(data(logical(data))),2,'significant');
    ymax = round(max(data(logical(data))),2,'significant');
else
    [ymin,ymax] = options.lims;
end
ymine = floor(log10(abs(ymin))); ymaxe = floor(log10(abs(ymax)));
ytlabelmin = num2str(ymin); ytlabelmax = num2str(ymax);
position(3) = position(3)/1.1;
if ymine
    ytlabelmin = [char(num2str(ymin/(10^ymine))) 'e' char(num2str(ymine))];
    % position(1) = position(1)+0.05;
    position(3) = position(3)*(1.15);
end
if ymaxe
    ytlabelmax = [char(num2str(ymax/(10^ymaxe))) 'e' char(num2str(ymaxe))];
    % position(1) = position(1)+0.05;
    position(3) = position(3)*(1.15);
end
ytlabel = {['\color[rgb]{1 1 1}' ytlabelmin],['\color[rgb]{1 1 1}' ytlabelmax]};

%
y = linspace(ymin, ymax, 100)'; xmin = 0; xmax = 1;
left = xmin*ones(size(y)); right = xmax*ones(size(y));
% Define the vertices: the points at (x, f(x)) and (x, 0)
N = length(y);
%verts = [x(:), top(:); x(:) bottom(:)];
verts = [left(:), y(:); right(:), y(:)];
% Define the faces to connect each adjacent f(x) and the corresponding points at y = 0.
q = (1:N-1)';
faces = [q, q+1, q+N+1, q+N];
% colormap
if ymin >= 0
    color1 = linspace(0.05*colormax(1),colormax(1),N)';
    color2 = linspace(0.05*colormax(2),colormax(2),N)';
    color3 = linspace(0.05*colormax(3),colormax(3),N)';
else
    % position(3) = position(3)*1.1;
    color1a = linspace(colormin(1),0,N/2)'; color1b = linspace(0,colormax(1),N/2)'; color1 = [color1a; color1b];
    color2a = linspace(colormin(2),0,N/2)'; color2b = linspace(0,colormax(2),N/2)'; color2 = [color2a; color2b];
    color3a = linspace(colormin(3),0,N/2)'; color3b = linspace(0,colormax(3),N/2)'; color3 = [color3a; color3b];
end

M = [color1, color2, color3];
ax = axes();
patch('Faces', faces, 'Vertices', verts, 'FaceVertexCData', [M; M], ...
    'FaceColor', 'interp', 'EdgeColor', 'none')
set(ax,'OuterPosition',position,...
        'YLim',[ymin ymax],'YTick',[ymin ymax],'XTick',[],...
        'YTickLabel',ytlabel,...
        'Color','k','xcolor','k','ycolor','k',...
        'XColor','k','YColor','k','FontSize',fontsize);
    
end