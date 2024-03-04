classdef circularGraph_lsd < handle
% CIRCULARGRAPH Plot an interactive circular graph to illustrate connections in a network.
%
%% Syntax
% circularGraph(X)
% circularGraph(X,'PropertyName',propertyvalue,...)
% h = circularGraph(...)
%
%% Description
% A 'circular graph' is a visualization of a network of nodes and their
% connections. The nodes are laid out along a circle, and the connections
% are drawn within the circle. Click on a node to make the connections that
% emanate from it more visible or less visible. Click on the 'Show All'
% button to make all nodes and their connections visible. Click on the
% 'Hide All' button to make all nodes and their connections less visible.
%
% Required input arguments.
% X : A symmetric matrix of numeric or logical values.
%
% Optional properties.
% Colormap : A N by 3 matrix of [r g b] triples, where N is the 
%            length(adjacenyMatrix).
% Label    : A cell array of N strings.
%%
% Copyright 2014 The MathWorks, Inc.
  properties
    Node = node(0,0); % Array of nodes
    ColorMap;         % Colormap for nodes
    Label;            % Cell array of strings
    Areacolors;       % Colors for dots
    SelfConnections;  % diagonal values
    Limits;           % sets limits on alpha and linewidths
    FontSize;         % sets the node label font size
    ShowButton;       % Turn all nodes on
    HideButton;       % Turn all nodes off
  end
  
  methods
    function this = circularGraph_lsd(adjacencyMatrix,varargin)
      % Constructor
      p = inputParser;
      
      defaultColorMap = repmat(255*[0.5 0.5 0.5],[length(adjacencyMatrix) 1]);
      defaultLabel = cell(length(adjacencyMatrix));
      defaultAreacolors = repmat([1 1 1],[length(adjacencyMatrix) 1]);
      defaultLimits = [min(reshape(adjacencyMatrix,[1,length(adjacencyMatrix)^2])),...
                       max(reshape(adjacencyMatrix,[1,length(adjacencyMatrix)^2]))];
      defaultFontSize = 14;
      for i = 1:length(defaultLabel)
        defaultLabel{i} = num2str(i);
      end
      defaultSelfConnections = ones([1,length(adjacencyMatrix)]);
      
      addRequired(p,'adjacencyMatrix',@(x)(isnumeric(x) || islogical(x)));
      addParameter(p,'ColorMap',defaultColorMap,@(colormap)length(colormap) == length(adjacencyMatrix));
      addParameter(p,'Label'   ,defaultLabel   ,@iscell);
      addParameter(p,'Areacolors',defaultAreacolors,@(colormap)length(colormap) == length(adjacencyMatrix));
      addParameter(p,'SelfConnections',defaultSelfConnections,@(x)(isnumeric(x)));
      addParameter(p,'Limits',defaultLimits,@(x)(isnumeric(x)));
      addParameter(p,'FontSize',defaultFontSize,@(x)(isnumeric(x)));
      
      parse(p,adjacencyMatrix,varargin{:});
      this.ColorMap = p.Results.ColorMap;
      this.Label    = p.Results.Label;
      this.Areacolors = p.Results.Areacolors;
      this.SelfConnections = p.Results.SelfConnections;
      this.Limits = p.Results.Limits;
      this.FontSize = p.Results.FontSize;
      
      %{
      this.ShowButton = uicontrol(...
        'Style','pushbutton',...
        'Position',[0 40 80 40],...
        'String','Show All',...
        'Callback',@circularGraph.showNodes,...
        'UserData',this);
      
      this.HideButton = uicontrol(...
        'Style','pushbutton',...
        'Position',[0 0 80 40],...
        'String','Hide All',...
        'Callback',@circularGraph.hideNodes,...
        'UserData',this);
      %}
      
      fig = gcf;
      set(fig,...
        'UserData',this,...
        'CloseRequestFcn',@circularGraph_lsd.CloseRequestFcn);
      
      % Draw the nodes
      delete(this.Node);
      t = linspace(-pi,pi,length(adjacencyMatrix) + 1).'; % theta for each node
      extent = zeros(length(adjacencyMatrix),1);
      
      % Node size based on self-connection value
      normselfconn = this.SelfConnections;
      normselfconn(~this.SelfConnections) = 0.4;
      normselfconn(logical(this.SelfConnections)) = ...
        normalize(-this.SelfConnections(logical(this.SelfConnections)),...
        'range',[0.5,1]);
      
      for i = 1:length(adjacencyMatrix)
        this.Node(i) = node(cos(t(i)),sin(t(i)));
        % this.Node(i).Color = this.ColorMap(i,:);
        this.Node(i).Label = this.Label{i};
        this.Node(i).NodeMarker.Color = this.Areacolors(i,:);
        % this.Node(i).TextLabel.FontSize = this.FontSize;

        %{
        % If increase (i.e. selfconn > 0)
        if this.SelfConnections(i) > 0
            % this.Node(i).NodeMarker.Marker = '^';
            this.Node(i).NodeMarker.MarkerEdgeColor = ...
                [255 181 0]./255;
        elseif this.SelfConnections(i) < 0 
            % this.Node(i).NodeMarker.Marker = 'v';
            this.Node(i).NodeMarker.MarkerEdgeColor = ...
                [0 201 239]./255;
        else
            this.Node(i).NodeMarker.Marker = 'o';
        end
        
                
        % Variable MarkerSize
        markersize = get(gca,'Position');
        markersize = (9*markersize(3)/432)*normselfconn(i);
        % markersize = 9*markersize(3)/432;
        this.Node(i).NodeMarker.MarkerSize = markersize;
        %}
        
        % Set size of marker and label font based on axis size
        figunits    = get(gcf,'Units');
        axunits     = get(gca,'Units');
        
        figurepos   = get(gcf,'OuterPosition');
        minfigsize  = min([figurepos(3),figurepos(4)]);
        axespos     = get(gca,'OuterPosition');
        minaxsize   = min([axespos(3),axespos(4)]);
        
        if axunits == 'normalized'
            minaxsize = minaxsize*minfigsize;
        end
        
        labelsize  = 0.001*minaxsize;
        markersize = 20 * labelsize;
        fontsize   = 14 * labelsize;
        
        this.Node(i).NodeMarker.MarkerSize = markersize;
        this.Node(i).TextLabel.FontSize = fontsize;
      end
      
      % Find non-zero values of s and their indices
      [row,col,v] = find(adjacencyMatrix);
      
      % Calculate alpha values based on values of s (stored in v)
      al_min = 0.05;
      al_scalar = 0.95;
      al = abs(v)./max(abs(this.Limits));
      al = al_scalar*al + al_min;
      % al = normalize(abs(v),'range',[0.05,1]);
      vs = v';
      
      % Calculate line widths based on values of s (stored in v).
      lw_min  = 0.1;
      lw_scalar = 1.5;
      lw = abs(v)./max(abs(this.Limits));
      if sum(lw) == numel(lw) % all lines are the same width.
        lw = repmat(lw_min,numel(lw),1);
      else % lines of variable width.
        lw = lw_scalar*lw + lw_min;
      end
      
      %{
      % Draw connections on the Poincare hyperbolic disk.
      %
      % Equation of the circles on the disk:
      % x^2 + y^2 
      % + 2*(u(2)-v(2))/(u(1)*v(2)-u(2)*v(1))*x 
      % - 2*(u(1)-v(1))/(u(1)*v(2)-u(2)*v(1))*y + 1 = 0,
      % where u and v are points on the boundary.
      %
      % Standard form of equation of a circle
      % (x - x0)^2 + (y - y0)^2 = r^2
      %
      % Therefore we can identify
      % x0 = -(u(2)-v(2))/(u(1)*v(2)-u(2)*v(1));
      % y0 = (u(1)-v(1))/(u(1)*v(2)-u(2)*v(1));
      % r^2 = x0^2 + y0^2 - 1
      %}
      
      for j = 1:length(v)
        [~,lineorder] = sort(lw,'ascend');
        i = lineorder(j);
        
        % for non-self-connections (off-diagonals)
        if row(i) ~= col(i)
          
          % assign colour based on positive/negative
          if vs(i)>0
            color_i = [255 181 0]./255;           
          else
            color_i = [0 201 239]./255;            
          end
        
          if abs(row(i) - col(i)) - length(adjacencyMatrix)/2 == 0 
            % points are diametric, so draw a straight line
            u = [cos(t(row(i)));sin(t(row(i)))];
            v = [cos(t(col(i)));sin(t(col(i)))];
            Xi = [u(1);v(1)]; Yi = [u(2);v(2)];  
            
          else % points are not diametric, so draw an arc
            u  = [cos(t(row(i)));sin(t(row(i)))];
            v  = [cos(t(col(i)));sin(t(col(i)))];
            x0 = -(u(2)-v(2))/(u(1)*v(2)-u(2)*v(1));
            y0 =  (u(1)-v(1))/(u(1)*v(2)-u(2)*v(1));
            r  = sqrt(x0^2 + y0^2 - 1);
            thetaLim(1) = atan2(u(2)-y0,u(1)-x0);
            thetaLim(2) = atan2(v(2)-y0,v(1)-x0);
            
            if u(1) >= 0 && v(1) >= 0 
              % ensure the arc is within the unit disk
              theta = [linspace(max(thetaLim),pi,50),...
                       linspace(-pi,min(thetaLim),50)].';
            else
              theta = linspace(thetaLim(1),thetaLim(2)).';
            end
            
            Xi = r*cos(theta)+x0; Yi = r*sin(theta)+y0;
            
          end
          
          % Draw line
          this.Node(row(i)).Connection(end+1) = line(...
            Xi,...
            Yi,...
            'LineWidth', lw(i),...
            'Color', [color_i al(i)],...
            'PickableParts','none');
        end
      end
      
      axis image;
      ax = gca;
      for i = 1:length(adjacencyMatrix)
        extent(i) = this.Node(i).Extent;
      end
      extent = max(extent(:));
      ax.XLim = ax.XLim + extent*[-1 1];
      fudgeFactor = 1.75; % Not sure why this is necessary. Eyeballed it.
      ax.YLim = ax.YLim + fudgeFactor*extent*[-1 1];
      ax.Visible = 'off';
      ax.SortMethod = 'childorder';
      
      
      set(gca,'color','k')
      set(gca,'FontName','Helvetica','XColor','w','YColor','w')
      set(gcf,'inverthardcopy','off');
      set(gcf,'color','k');
    end
    
  end
  
  methods (Static = true)
    function showNodes(this,~)
      % Callback for 'Show All' button
      n = this.UserData.Node;
      for i = 1:length(n)
        n(i).Visible = true;
      end
    end
    
    function hideNodes(this,~)
      % Callback for 'Hide All' button
      n = this.UserData.Node;
      for i = 1:length(n)
        n(i).Visible = false;
      end
    end
    
    function CloseRequestFcn(this,~)
      % Callback for figure CloseRequestFcn
      c = this.UserData;
      for i = 1:length(c.Node)
        delete(c.Node(i));
      end
      delete(gcf);
    end
    
  end
  
end