function cognemo_ringlabel2(arealabels,options)
%% Preamble
%{
This function places a ring of region-area labels around the connectogram
plots
%}
%% Unpack options

radius = options.radius;
fontsize = options.fontsize;
currfig = gcf; boxsize = currfig.Position(3); % fontsize = fontsize*boxsize/8.6389;

%%

N_r = length(arealabels); d_arc = 2*pi/N_r;
% arealabels_unique = unique(arealabels,'stable');
% N_areas = length(arealabels_unique);

%% Initialize variable arrays for annotations

% N_occ_area = zeros([N_areas,1]);
% angle = zeros([N_areas,1]);
% H_al = repmat(" ",[N_areas,1]); V_al = repmat(" ",[N_areas,1]);
% X = zeros([N_areas,1]); Y = zeros([N_areas,1]);

%% Adjusting angles to match dot labels:
% the 'dot' labels are rotated back from [1,0] by pi
angle_adjust = pi/2; Q = zeros([4,2]);
Q(1,:) = [0,pi/2]; Q(2,:) = [pi/2,pi]; Q(3,:) = [pi,3*pi/2]; Q(4,:) = [3*pi/2,2*pi];
Q = Q + angle_adjust;

%% First iteration:
% Must set the first one to initiate recursive definition in loop
%{
N_occ_area = 1; j = 1;
for i = 1:N_r
    if arealabels(i+1)==arealabels(i)
        N_occ_area = N_occ_area + 1;
    else
        j = i+1;
        break
    end
end
angle = d_arc*floor(0.5*N_occ_area)+angle_adjust;
H_al = "Right"; V_al = "Top"; % Q1
X = radius*cos(angle); Y = radius*sin(angle);
% Draw annotation
text(X,Y,arealabels_unique,...
     'Color','w',...
     'FontSize',fontsize,...
     'FontSmoothing','on',...
     'HorizontalAlignment',H_al,...
     'VerticalAlignment',V_al);
%}
%% Cheat version
%{
i = 1; j = 1; N_occ_area = [];
while j < N_r
    
    % Set angles, alignment, and position
    N_occ_area_i = 1;
    for l = j+1:N_r
        if arealabels(l)==arealabels(l-1)
            N_occ_area_i = N_occ_area_i + 1;
            j = j+1;
        else
            N_occ_area = [N_occ_area N_occ_area_i];
            j = l+1;
            i = i+1;
            break
        end
    end
    angle_i = d_arc*floor(0.5*N_occ_area_i) + d_arc*sum(N_occ_area(1:i-1)) + angle_adjust;
    if (angle_i > Q(1,1)) && (angle_i < Q(1,2))
        H_al_i = "Right"; V_al_i = "Top"; % Q1
    elseif (angle_i > Q(2,1)) && (angle_i < Q(2,2))
        H_al_i = "Left"; V_al_i = "Top"; % Q2
    elseif (angle_i > Q(3,1)) && (angle_i < Q(3,2))
        H_al_i = "Left"; V_al_i = "Bottom"; % Q3
    elseif (angle_i > Q(4,1)) && (angle_i < Q(4,2))
        H_al_i = "Right"; V_al_i = "Bottom"; % Q4
    end
    X_i = radius*cos(angle_i); Y_i = radius*sin(angle_i);
        
    % Draw annotation
    text_i = text(X_i,Y_i,arealabels(j-1),...
                 'Color','w',...
                 'FontSize',fontsize,...
                 'FontSmoothing','on',...
                 'HorizontalAlignment',H_al_i,...
                 'VerticalAlignment',V_al_i);
     
     % add Rotation
     if isfield(options,'rotate')
         if options.rotate
             t = atan2(Y_i,X_i);
             
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
%}

arealabels = ["PFr",...
              "Fr",...
              "Ins",...
              "Tem",...
              "Par",...
              "Occ",...
              "SbC",...
              "CeB",...
              "Ver",...
              "CeB",...
              "SbC",...
              "Occ",...
              "Par",...
              "Tem",...
              "Ins",...
              "Fr"];

angle = [0,...
         0.21,...
         0.31,...
         0.4,...
         0.6,...
         0.78,...
         0.9,...
         0.95,...
         0.98,...
         1,...
         1.05,...
         1.16,...
         1.35,...
         1.54,...
         1.64,...
         1.76];
angle = pi*angle + angle_adjust;
    
for i = 1:length(arealabels)
    X_i = radius*cos(angle(i)); Y_i = radius*sin(angle(i));
    if angle(i) == angle_adjust
        H_al_i = "Center"; V_al_i = "Bottom";
    elseif (angle(i) > Q(1,1)) && (angle(i) < Q(1,2))
        H_al_i = "Right"; V_al_i = "Bottom"; % Q1
    elseif (angle(i) > Q(2,1)) && (angle(i) < Q(2,2))
        H_al_i = "Right"; V_al_i = "Top"; % Q2
    elseif (angle(i) > Q(3,1)) && (angle(i) < Q(3,2))
        H_al_i = "Left"; V_al_i = "Top"; % Q3
    elseif (angle(i) > Q(4,1)) && (angle(i) < Q(4,2))
        H_al_i = "Left"; V_al_i = "Bottom"; % Q4
    elseif angle(i) == pi + angle_adjust
        H_al_i = "Center"; V_al_i = "Top";
    end
    text(X_i,Y_i,arealabels(i),...
                     'Color','w',...
                     'FontSize',fontsize,...
                     'FontSmoothing','on',...
                     'HorizontalAlignment',H_al_i,...
                     'VerticalAlignment',V_al_i); clear X_i Y_i H_al_i V_al_i text_i
end
end