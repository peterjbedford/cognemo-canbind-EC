function cognemo_ringlabel(arealabels,options)
%% Preamble
%{
This function places a ring of region-area labels around the connectogram
plots
%}
%% Unpack options

radius = options.radius;
fontsize = options.fontsize;
currfig = gcf; boxsize = currfig.Position(3); %fontsize = fontsize*boxsize/8.6389;

%%

N_r = length(arealabels); d_arc = 2*pi/N_r;
arealabels_unique = unique(arealabels,'stable');
N_areas = length(arealabels_unique);

%% Initialize variable arrays for annotations

N_occ_area = zeros([N_areas,1]); angle = zeros([N_areas,1]);
H_al = repmat(" ",[N_areas,1]); V_al = repmat(" ",[N_areas,1]);
X = zeros([N_areas,1]); Y = zeros([N_areas,1]);

%% Adjusting angles to match dot labels:
% the 'dot' labels are rotated back from [1,0] by pi
angle_adjust = -pi; Q = zeros([4,2]);
Q(1,:) = [0,pi/2]; Q(2,:) = [pi/2,pi]; Q(3,:) = [pi,3*pi/2]; Q(4,:) = [3*pi/2,2*pi];
Q = Q + angle_adjust;

%% First iteration:
% Must set the first one to initiate recursive definition in loop

N_occ_area(1) = length(find(arealabels==arealabels_unique(1)));
angle(1) = d_arc*floor(0.5*N_occ_area(1))+angle_adjust;
H_al(1) = "Right"; V_al(1) = "Top"; % Q1
X(1) = radius*cos(angle(1)); Y(1) = radius*sin(angle(1));
% Draw annotation
text(X(1),Y(1),arealabels_unique(1),...
     'Color','w',...
     'FontSize',fontsize,...
     'FontSmoothing','on',...
     'HorizontalAlignment',H_al(1),...
     'VerticalAlignment',V_al(1));
 
%% Loop
for i = 2:N_areas
    
    % Set angles, alignment, and position
    N_occ_area(i) = length(find(arealabels==arealabels_unique(i)));
    angle(i) = d_arc*floor(0.5*N_occ_area(i)) + d_arc*sum(N_occ_area(1:i-1)) + angle_adjust;
    if (angle(i) > Q(1,1)) && (angle(i) < Q(1,2))
        H_al(i) = "Right"; V_al(i) = "Top"; % Q1
    elseif (angle(i) > Q(2,1)) && (angle(i) < Q(2,2))
        H_al(i) = "Left"; V_al(i) = "Top"; % Q2
    elseif (angle(i) > Q(3,1)) && (angle(i) < Q(3,2))
        H_al(i) = "Left"; V_al(i) = "Bottom"; % Q3
    elseif (angle(i) > Q(4,1)) && (angle(i) < Q(4,2))
        H_al(i) = "Right"; V_al(i) = "Bottom"; % Q4
    end
    X(i) = radius*cos(angle(i)); Y(i) = radius*sin(angle(i));
        
    % Draw annotation
    text_i = text(X(i),Y(i),arealabels_unique(i),...
                 'Color','w',...
                 'FontSize',fontsize,...
                 'FontSmoothing','on',...
                 'HorizontalAlignment',H_al(i),...
                 'VerticalAlignment',V_al(i));
     
     % add Rotation
     if isfield(options,'rotate')
         if options.rotate
             t = atan2(Y(i),X(i));
             
             text_i.VerticalAlignment = 'middle';
             
             if abs(t) > pi/2
                text_i.Rotation = 180*(t/pi + 1);
                text_i.HorizontalAlignment = 'right';
                
             else
                text_i.Rotation = t*180/pi;
                text_i.HorizontalAlignment = 'left';
             end
         end
     end
     
end; clear i

end